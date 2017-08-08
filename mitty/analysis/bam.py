import time
import logging

import pysam

from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname

logger = logging.getLogger(__name__)


def is_single_end_bam(bam_fname):
  bam_fp = pysam.AlignmentFile(bam_fname)
  r = next(bam_fp, None)
  return not r.is_paired if r is not None else True  # Empty BAM? Don't care


def bam_iter(bam_fname, sidecar_fname, limit=None):
  """Given a BAM file path return us tuples of paired reads.

  :param bam_fname: BAM file name
  :param sidecar_fname: long qname overflow file
  :param limit: If not None limit the number of reads scanned to this count
  :return: a generator that returns pairs of reads from the file
  """
  se_bam = is_single_end_bam(bam_fname)
  long_qname_table = load_qname_sidecar(sidecar_fname)

  read_dict = {}
  for n, rd in enumerate( pysam.AlignmentFile( bam_fname ).fetch( until_eof=True )):
    if limit is not None and n >= limit:
      break
    if rd.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    ri = parse_qname(rd.qname, long_qname_table=long_qname_table)[1 if rd.is_read2 else 0]
    if se_bam:
      yield (rd, ri, True)
    else:
      key = rd.qname[:20]
      if key not in read_dict:
        read_dict[key] = [None, None]

      rl = read_dict[key]
      rl[0 if rd.is_read1 else 1] = (rd, ri, True)

      if all(rl):
        yield rl
        del read_dict[key]


def derr(r_iter, d_max):
  """Return reads with XD tag filled out with d_err metric

  :param r_iter: An iterable of read tuples
  :param d_max:
  :return:
  """
  for r in r_iter:
    for mate in r:
      d_err = score_alignment_error(r=mate[0], ri=mate[1], max_d=d_max)
      mate[0].set_tag('XD', d_err)
    yield r


def accept_reads(r_iter, f):
  for r in r_iter:
    new_r = [
      (mate[0], mate[1], f(mate) and mate[2])
      for mate in r
    ]
    if any([mate[2] for mate in new_r]):
      yield new_r


def discard_ref(r_iter):
  """Discard reference reads

  :param r_iter:
  :return:
  """
  for r in accept_reads(r_iter, lambda mate: len(mate[1].v_list) > 0):
    yield r


def discard_non_ref(r_iter):
  """Discard non-reference reads

  :param r_iter:
  :return:
  """
  for r in accept_reads(r_iter, lambda mate: len(mate[1].v_list) == 0):
    yield r


def discard_derr(r_iter, d_range):
  """Discard reads falling within given d_range

  :param r_iter:
  :param d_range: (low_d_err, high_d_err) e.g. (-1000, -10) or (10, 1000)
  :return:
  """
  for r in accept_reads(r_iter, lambda mate: not (d_range[0] <= mate[0].get_tag('XD') <= d_range[1])):
    yield r


def discard_v(r_iter, v_range):
  """Discard reads with variants falling within given v_range

  :param r_iter:
  :param v_range: (low_v_size, high_v_size) e.g. (-1000, -50) or (50, 1000) or (-50, 50)
  :return:
  """
  for r in accept_reads(r_iter,
                        lambda mate: not all(v_range[0] <= v <= v_range[1]
                                         for v in mate[1].v_list)):
    yield r

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
