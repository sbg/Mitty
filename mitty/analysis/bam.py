import time
import logging
from collections import namedtuple, Counter

import numpy as np
import pysam

from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname

logger = logging.getLogger(__name__)

Read = namedtuple('Read',
                  ['read', 'read_info', 'd_err', 'filter_pass', 'cat_list'])


def is_single_end_bam(bam_fname):
  bam_fp = pysam.AlignmentFile(bam_fname)
  r = next(bam_fp, None)
  return not r.is_paired if r is not None else True  # Empty BAM? Don't care


# Data Sources ----------------------------------------------------------------


def bam_iter(bam_fname, sidecar_fname, max_d=200, every=None):
  """Given a BAM file path return us tuples of paired reads.

  :param bam_fname: BAM file name
  :param sidecar_fname: long qname overflow file
  :param max_d: maximum d_err we consider
  :param every: If not None, returns every Nth read/read-pair
  :return: a generator that returns pairs of reads from the file
  """
  se_bam = is_single_end_bam(bam_fname)
  long_qname_table = load_qname_sidecar(sidecar_fname)

  read_dict = {}
  ctr = 0
  for rd in pysam.AlignmentFile(bam_fname).fetch(until_eof=True):
    if rd.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    ri = parse_qname(rd.qname, long_qname_table=long_qname_table)[1 if rd.is_read2 else 0]
    d_err = score_alignment_error(r=rd, ri=ri, max_d=max_d)

    if se_bam:
      key = None
      rl = [None]
    else:
      key = rd.qname[:20]
      if key not in read_dict:
        read_dict[key] = [None, None]
      rl = read_dict[key]

    rl[0 if (rd.is_read1 or se_bam) else 1] = Read(rd, ri, d_err, True, [])

    if all(rl):
      if every is None or ctr == 0:
        yield rl
        ctr = every or 0
      if key is not None:
        del read_dict[key]
      ctr -= 1


# Data sinks ------------------------------------------------------------------


def simple_sink(r_iter):
  """This just consumes the reads so that our filter chain is processed. Comes in useful in
  some cases when we just drop the reads off at the end of the pipeline

  :param r_iter:
  :return:
  """
  for r in r_iter:
    pass


def write_to_bam(r_iter, fname, header):
  out_fp = pysam.AlignmentFile(fname, 'wb', header=header)
  for r in r_iter:
    for mate in r:
      out_fp.write(mate.read)


# Filtering operations --------------------------------------------------------


def tag_derr(r_iter):
  """Return reads with XD tag filled out with d_err metric

  :param r_iter: An iterable of read tuples
  :return:
  """
  for r in r_iter:
    for mate in r:
      mate.read.set_tag('XD', mate.d_err)
    yield r


def filter_reads(r_iter, f):
  """Filter out reads based on f

  :param r_iter:
  :param f: filter function
  :return:
  """
  for r in r_iter:
    new_r = [
      mate._replace(filter_pass=f(mate) and mate.filter_pass)
      for mate in r
    ]
    if any([mate.filter_pass for mate in new_r]):
      yield new_r


def categorize_reads(r_iter, f_dict):
  """Fill in cat_list of Reads. Note that there is no loss of reads in this function.
  If a read does not match any filter cat_list is empty, which corresonds to 'uncategorized'

  :param r_iter:
  :param f_dict: Dictionary of key: filter_function pairs
                 key will go into cat_list field of read if filter passes
  :return:
  """
  for r in r_iter:
    yield [
      mate._replace(cat_list=[k for k, f in f_dict.items() if f(mate)])
      for mate in r
    ]


def count_reads(r_iter, counter=None):
  """

  :param r_iter:
  :param counter: a dictionary of counts, can be empty
  :return:
  """
  for r in r_iter:
    for mate in r:
      for cat in (mate.cat_list or ['nocat']):
        if cat not in counter:
          counter[cat] = 0
        counter[cat] += 1
    yield r


def count_reads_sink(r_iter):
  """This is a sink. It does not act as a filter

  :param r_iter:
  :return:
  """
  return Counter(
    cat
    for r in r_iter
    for mate in r
    for cat in (mate.cat_list or ['nocat'])
    if mate.filter_pass
  )


def zero_dmv(max_d=200, max_MQ=70, max_vlen=200):
  return np.zeros(shape=(2 * max_d + 3, max_MQ + 1, 2 * max_vlen + 1 + 2 + 1), dtype=int)


def alignment_hist(r_iter, dmv_mat):
  """Compute the dmv matrix which is as defined as follows:

  A 3D matrix with dimensions:
    Xd - alignment error  [0]  -max_d, ... 0, ... +max_d, wrong_chrom, unmapped (2 * max_d + 3)
    MQ - mapping quality  [1]  0, ... max_MQ (max_MQ + 1)
    vlen - length of variant carried by read [2]  -max_vlen, ... 0, ... +max_vlen, Ref, Margin
                                                  (2 * max_vlen + 1 + 2)

  * The ends of the ranges (-max_d, max_d) (-max_vlen, +max_vlen) include that value and all values exceeding
  * Ref collects all the reference reads
  * Margin collects the marginal sums for d x MQ. This is necessary because one read can appear in multiple
    bins on the vlen axis and cause a slight excess of counts when marginals are computed

  :param r_iter:
  :param dmv_mat:
  :return: iterator
  """
  max_d = int((dmv_mat.shape[0] - 3) / 2)
  max_MQ = int(dmv_mat.shape[1] - 1)
  max_vlen = int((dmv_mat.shape[2] - 4) / 2)

  for r in r_iter:
    for mate in r:
      i = max_d + mate.d_err
      j = min(mate.read.mapping_quality, max_MQ)
      if mate.read_info.v_list:
        for v_size in mate.read_info.v_list:
          k = min(max_vlen, max(0, max_vlen + v_size))
          dmv_mat[i, j, k] += 1
      else:
        k = 2 * max_vlen + 1
        dmv_mat[i, j, k] += 1

    yield r


# Library of useful filter functions ------------------------------------------


def non_ref():
  """Keep reads with variants

  :param mate: Read object
  :return: T/F
  """
  return lambda mate: len(mate.read_info.v_list) > 0


def pure_ref():
  """Keep reads with no variants

  :param mate: Read object
  :return: T/F
  """
  return lambda mate: len(mate.read_info.v_list) == 0


def discard_derr(d_range):
  """Discard reads falling within given d_range

  :param d_range: (low_d_err, high_d_err) e.g. (-1000, -10) or (10, 1000)
  :return:
  """
  return lambda mate: not (d_range[0] <= mate.d_err <= d_range[1])


def discard_v(v_range):
  """Discard reads with variants falling within given v_range

  :param v_range: (low_v_size, high_v_size) e.g. (-1000, -50) or (50, 1000) or (-50, 50)
  :return:
  """
  return lambda mate: not all(v_range[0] <= v <= v_range[1]
                              for v in mate.read_info.v_list)




"""

`derr ( r, d_max )` - return reads with XD tag filled out with d_err metric
`discard_ref ( r )` - discard reference reads
`discard_non_ref ( r )` - discard non-reference reads
`filter_derr ( r, d_range, d_max )` - given a read iterator filter out reads falling in given d_range
`filter_v ( r, v_range )` - filter out reads with no variants outside this range
`categorize ( r, cat_dict )` - given a dictionary of filter functions create categories (or sub-categories)
`count ( r, result)` - count up all the reads for each category, return in the result dictionary
`save_to_bam( r, header )` - save reads in separate BAM files with names matching the category names
`read_fate_plot( result, ax )` - plot a histogram with different category labels based on a result dictionary
                                 passing multiple results will result in multiple histograms on the same axes
`xmv ( r, d_max, MQ_max, vlen_max, result )` - Three dimensional alignment analysis histogram bin counts
`xmv_plot ( result, ax )` - alignment analysis plot. Multiple plots on same axes if

"""
