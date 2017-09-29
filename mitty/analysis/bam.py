import pickle
import gzip
import time
from collections import OrderedDict
import logging

import cytoolz
import numpy as np
import pysam
# import pandas as pd
import xarray as xr

from mitty.benchmarking.alignmentscore import score_alignment_error, correct_tlen, load_qname_sidecar, parse_qname

logger = logging.getLogger(__name__)

# Using a dict which is a builtin and pretty fast
# https://gist.github.com/cordella/2861038
# read_dict = {
#   'read': ...
#   'read_info': ...
#   'd_err': ...
#   'cat_list': ...
#   'fpass': ...
# }


def get_header(bam_fname):
  return pysam.AlignmentFile(bam_fname).header


# Data Sources ----------------------------------------------------------------


def read_bam(bam_fname):
  """Fetch reads (including unmapped) from a BAM.

  :param bam_fname: path to BAM file
  :param sidecar_fname: If provided, will fill out read_info
  :param every: Only take every nth read
  :return:
  """
  for rd in pysam.AlignmentFile(bam_fname).fetch(until_eof=True):
    if rd.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    yield (
      {
        'read': rd,
        'fpass': True
      },
    )


# Data sinks ------------------------------------------------------------------


def simple_sink(riter):
  """This just consumes the reads so that our filter chain is processed. Comes in useful in
  some cases when we just drop the reads off at the end of the pipeline

  :param riter:
  :return:
  """
  for r in riter:
    pass


@cytoolz.curry
def write_bam(fname, header, riter):
  out_fp = pysam.AlignmentFile(fname, 'wb', header=header)
  for r in riter:
    for mate in r:
      out_fp.write(mate['read'])
    yield r


# Filtering operations --------------------------------------------------------


@cytoolz.curry
def filter_reads(f, condition, riter):
  """Filter out reads based on f

  :param f: filter function
  :param condition: is either one of the python functions all or any
                    to indicate if we should accept paired reads only if both the mates pass
                    or if any of the mates pass
  :param riter:
  :return:
  """
  for r in riter:
    # TODO: looks like we don't need 'fpass'
    new_r = tuple(dict(mate, fpass=f(mate) and mate['fpass']) for mate in r)
    if condition(tuple(mate['fpass'] for mate in new_r)):
      yield new_r


# Library of useful filter functions ------------------------------------------


def non_ref():
  """Keep reads with variants

  :return: function that takes mate as input and returns T/F
  """
  return lambda mate: len(mate['read_info'].v_list) > 0


def pure_ref():
  """Keep reads with no variants

  :return: function that takes mate as input and returns T/F
  """
  return lambda mate: len(mate['read_info'].v_list) == 0


def derr(min, max):
  """Keep reads within this limit

  :param min:
  :param max:
  :return:
  """
  return lambda mate: min <= mate['d_err'] <= max


def vsize(min, max):
  """Keep reads with at least one variant falling within given v_range
  Returns False for pure reference reads

  :param min:
  :param max:
  :return:
  """
  return lambda mate: any(min <= v <= max for v in mate['read_info'].v_list)


# Processing tools ------------------------------------------------------------


def make_pairs(riter):
  """Given a stream of reads pair them up by qname. Assumes riter is spewing tuples of
  size 1

  :param riter:
  :return:
  """
  read_dict = {}
  for r in riter:
    rd = r[0]['read']
    key = rd.qname[:20]
    if key not in read_dict:
      read_dict[key] = r[0]
    else:  # You complete me
      yield (r[0], read_dict[key]) if rd.is_read1 else (read_dict[key], r[0])
      del read_dict[key]


@cytoolz.curry
def parse_read_qnames(sidecar_fname, titer):
  """Mutates dictionary: adds 'read_info' field to it.

  :param titer:
  :return:
  """
  long_qname_table = load_qname_sidecar(sidecar_fname) if sidecar_fname is not None else None

  for tplt in titer:
    ri = parse_qname(
        tplt[0]['read'].qname,
        long_qname_table=long_qname_table
    ) if long_qname_table is not None else [None, None]
    for mate in tplt:
      rd = mate['read']
      mate['read_info'] = ri[1 if rd.is_read2 else 0]
    yield tplt


@cytoolz.curry
def compute_derr(riter, max_d=200):
  """Mutates dictionary: adds d_err field to it.

  :param max_d:
  :param riter:
  :return:
  """
  for r in riter:
    for mate in r:
      mate['d_err'] = score_alignment_error(r=mate['read'], ri=mate['read_info'], max_d=max_d)
    yield r


@cytoolz.curry
def categorize_reads(f_dict, r_iter):
  """Fill in cat_list of Reads. Note that there is no loss of reads in this function.
  If a read does not match any filter cat_list is empty, which corresponds to 'uncategorized'
  Mutates read dictionary: adds 'cat_list' field to it. If a 'cat_list' field already exists
  it appends to it.

  :param r_iter:
  :param f_dict: Dictionary of key: filter_function pairs
                 key will go into cat_list field of read if filter passes

  e.g. f_dict = {
    'd = 0': lambda mate: mate['d_err'] == 0,
    '0 < d <= 50': lambda mate: 0 < mate['d_err'] <= 50,
    '50 < d': lambda mate: 50 < mate['d_err'] < 200,
    'WC': lambda mate: 200 < mate['d_err'],
    'UM': lambda mate: mate['read'].is_unmapped
    }

  :return: iterator
  """
  for r in r_iter:
    for mate in r:
      mate['cat_list'] = mate.get('cat_list', []) + [k for k, f in f_dict.items() if f(mate)]
    yield r


@cytoolz.curry
def count_reads(counter, r_iter):
  """The reads need to have gone through `categorize_reads` so that they have the
    'cat_list'] field filled out

  :param counter: a dictionary of counts
  :param r_iter:
  :return:
  """
  for r in r_iter:
    for mate in r:
      for cat in (mate['cat_list'] or ['nocat']):
        if cat not in counter:
          counter[cat] = 0
        counter[cat] += 1
    yield r


def fastmultihist(sample, bins):
  """

  :param sample:
  :param bins:
  :return:

  Because np.histogramdd is broken in terms of performance.
  Because we don't need no nanny state
  """
  N, D = sample.shape

  # Compute the bin number each sample falls into.
  Ncount = [np.digitize(sample[:, i], bins[i]) - 1 for i in np.arange(D)]

  # Unravel the D-dimensional coordinates onto a 1-D index
  nbin = np.array([len(b) - 1 for b in bins], int)
  xy = np.zeros(N, int)
  for i in np.arange(0, D - 1):
    xy += Ncount[i] * nbin[i + 1:].prod()
  xy += Ncount[-1]

  # Histogram the indices and then ravel the index counts back into D-dimensions
  return np.bincount(xy, minlength=nbin.prod()).reshape(nbin)


def make_tick_labels(bin_edges, var_label='x'):
  b = bin_edges
  return [
    '{} <= {} < {}'.format(round(b[i], 2), var_label, round(b[i+1], 2))
    for i in range(len(b) - 1)
  ]


def initialize_pah(
  xd_bin_edges=(-np.inf, -100, -40, -20, -10, -5, 0, 1, 6, 11, 21, 41, 101, np.inf),
  max_d=200,
  mq_bin_edges=(0, 1, 11, 41, 51, np.inf),
  v_bin_edges=(-np.inf, -20, 0, 1, 21, np.inf),
  xt_bin_edges=(-np.inf, -100, -40, -20, -10, -5, 0, 1, 6, 11, 21, 41, 101, np.inf),
  t_bin_edges=(-np.inf, 200, 301, 401, 601, np.inf),
  name='PairedAlignmentHistogram'):
  """Hijacks xarray.DataArray to store the histogram + some metadata

  The 'coords' field stores the actual spatial values (bin centers) you may want to use for plotting
  while the attrs dictionary stores the bin labels (which reminds you of any special bins there may exist)
  and descritpions of the

  :param xd_bin_edges:
  :param max_d:
  :param mq_bin_edges:
  :param v_bin_edges:
  :param xt_bin_edges: template length error bins
  :param t_bin_edges: correct template length bins
  :param buf_size: compute the histogram in chunks of N

  Binning rule is x0 <= x < x1

  The size of the matrix using the defaults is

  15 * 15 * 8 * 8 * 10 * 10 * 13 * 8 =  150 million elements

  """
  # Some of the dimensions need extra bins for particular exceptions to the data
  # xd_bin_edges has to be adjusted to accommodate WC and UM reads
  xd_bin_edges = list(xd_bin_edges)
  xd_bin_edges[-1] = max_d + 1
  xd_bin_edges += [max_d + 2, max_d + 3]
  # ...., max_d + 1, max_d + 2, max_d + 3
  #     |          |          \----- UM
  #     |          \----- WC
  #     \---- max_d+
  #
  # Also note that alignment analysis (not us) takes care of clipping d_err at max_d
  xd_labels = make_tick_labels(xd_bin_edges, 'xd')
  xd_labels[-2] = 'WC'
  xd_labels[-1] = 'UM'

  max_v = v_bin_edges[-2] + 1  # This is what we clip v size to be
  # xd_bin_edges has to be adjusted to accommodate Ref and V+
  v_bin_edges = list(v_bin_edges)
  v_bin_edges[-1] = max_v + 1
  v_bin_edges += [max_v + 2, max_v + 3]
  # ...., max_v + 1, max_v + 2, max_v + 3
  #     |          |          \-- V+
  #     |          \-- Ref
  #     \-- varlen+
  v_labels = make_tick_labels(v_bin_edges, 'V')
  v_labels[-2] = 'Ref'
  v_labels[-1] = 'V+'

  dimensions = [
      ('xd1', 'alignment error mate 1', np.array(xd_bin_edges), xd_labels),
      ('xd2', 'alignment error mate 2', np.array(xd_bin_edges), xd_labels),
      ('mq1', 'mapping quality mate 1', np.array(mq_bin_edges), make_tick_labels(mq_bin_edges, 'MQ')),
      ('mq2', 'mapping quality mate 2', np.array(mq_bin_edges), make_tick_labels(mq_bin_edges, 'MQ')),
      ('v1', 'length of variant carried by mate 1', np.array(v_bin_edges), v_labels),
      ('v2', 'length of variant carried by mate 2', np.array(v_bin_edges), v_labels),
      ('xt', 'template length error', np.array(xt_bin_edges), make_tick_labels(xt_bin_edges, 'Xt')),
      ('t', 'Correct template length', np.array(t_bin_edges), make_tick_labels(t_bin_edges, 'T'))
  ]

  dimension_metadata = OrderedDict(
    [
      (k[0],
       {
         'description': k[1],
         'bin_edges': k[2],
         'bin_centers': (k[2][:-1] + k[2][1:]) / 2
       })
      for k in dimensions
    ]
  )

  metadata = {
    'max_d': max_d,
    'max_v': max_v,
    'description': [k[1] for k in dimensions],
    'bin_centers': [(k[2][:-1] + k[2][1:]) / 2 for k in dimensions],
    'bin_edges': [k[2] for k in dimensions]  # This is bin edges in the form of a list so that histogram can use it
  }

  return xr.DataArray(
    data=np.zeros(
      shape=tuple(len(d[3]) for d in dimensions),
      dtype=int),
    coords=OrderedDict([(d[0], d[3]) for d in dimensions]),
    dims=[d[0] for d in dimensions],
    name=name,
    attrs=metadata
  )


@cytoolz.curry
def histogramize(titer, pah=None, buf_size=1000000):
  """Iterate over reads and fill out the histogram.

  This works as a sink. For some reason, if I write this as a filter
  nothing after the yield is executed and so we can't do the final finalize :/

  :return:
  """
  def _parse_vlist(_vl):
    if _vl:
      if len(_vl) == 1:
        return max(-max_v, min(_vl[0], max_v))  # Single variant read
      else:
        return max_v + 2  # Multi variant read
    else:
      return max_v + 1  # Reference read

  def _finalize():
    """given a buffer of data bin it into the histgram
    """
    pah.data += fastmultihist(buf[:idx, :], pah.attrs['bin_edges'])

  if pah is None:
    pah = initialize_pah()

  max_d, max_v = pah.attrs['max_d'], pah.attrs['max_v']
  temp_read = pysam.AlignedSegment()

  bs = buf_size
  buf = np.empty(shape=(bs, len(pah.coords)))
  idx = 0

  for tpl in titer:
    # For efficiency this assumes we have PE reads and have paired them up using make_pairs
    ctl = correct_tlen(tpl[0]['read_info'], tpl[1]['read_info'], temp_read)

    # Xd1 x Xd2 x MQ1 x MQ2 x vl1 x vl2 x Xt
    buf[idx, :] = [
      tpl[0]['d_err'],
      tpl[1]['d_err'],
      tpl[0]['read'].mapping_quality,
      tpl[1]['read'].mapping_quality,
      _parse_vlist(tpl[0]['read_info'].v_list),
      _parse_vlist(tpl[1]['read_info'].v_list),
      abs(tpl[0]['read'].template_length) - ctl,
      ctl
    ]
    idx += 1

    if idx >= bs:
      _finalize()
      idx = 0  # We don't bother to clear buf - we just over write it

  _finalize()


def save_pah(pah, fname):
  """Saves gzipped to pickle. Nice compression, imperceptibly slower than unzipped version but
  a LOT smaller. xarray's native serialization formats are hillariously broken.

  :param pah:
  :param fname:
  :return:
  """
  # from a hint http://henrysmac.org/blog/2010/3/15/python-pickle-example-including-gzip-for-compression.html
  with gzip.open(fname, 'wb') as f:
    pickle.dump(pah, f, protocol=-1)


def load_pah(fname):
  with gzip.open(fname, 'rb') as f:
    return pickle.load(f)


# TODO: copy over the attrs as well
def collapse(pah, **m):
  """Returns the marginal for the dimensions supplied.

  First we slice the matrix where we want
  Then we marginalize over the dimensions we want

  :param pah: the histogram
  :param m: each keyword argument must match a dimension
            if the value of the argument is None, then we keep that dimension
            if the value is a tuple (i, j), we slice out i:j from that dimension and then marginalize
            dimensions not mentioned are fully marginalized
  :return:
  """
  assert all(n in pah.dims for n in m.keys()), 'Mistake in dimension names'
  slices = tuple(slice(*(m.get(k, None) or [None])) for k in pah.dims)
  # slices = [..., slice(None), slice(i, j), .... ]
  # here slice(None) is for both the margin dimensions as well as the fully marginalized dimensions
  # slice(i, j) is for the partially marginalized dimensions
  collapse_idx_l = tuple(i for i in range(len(pah.dims))
                         if pah.dims[i] not in [d for d, v in m.items() if v is None])
  return pah[slices].sum(axis=collapse_idx_l)


def zero_dmv(max_d=200, max_MQ=70, max_vlen=200):
  return np.zeros(shape=(2 * max_d + 1 + 2, max_MQ + 1, 2 * max_vlen + 1 + 2), dtype=int)


@cytoolz.curry
def alignment_hist(dmv_mat, riter):
  """Compute the dmv matrix which is as defined as follows:

  A 3D matrix with dimensions:
    Xd - alignment error  [0]  -max_d, ... 0, ... +max_d, wrong_chrom, unmapped (2 * max_d + 1 + 2)
    MQ - mapping quality  [1]  0, ... max_MQ (max_MQ + 1)
    vlen - length of variant carried by read [2]  -max_vlen, ... 0, ... +max_vlen, Ref, Margin, Multiple
                                                  (2 * max_vlen + 1 + 2)



  * The ends of the ranges (-max_d, max_d) (-max_vlen, +max_vlen) include that value and all values exceeding
  * Ref collects all the reference reads
  * Margin collects the marginal sums for d x MQ. This is necessary because one read can appear in multiple
    bins on the vlen axis and cause a slight excess of counts when marginals are computed

  :param dmv_mat:
  :param riter:
  :return: iterator
  """
  max_d = int((dmv_mat.shape[0] - 3) / 2)
  max_MQ = int(dmv_mat.shape[1] - 1)
  max_vlen = int((dmv_mat.shape[2] - 3) / 2)

  for r in riter:
    for mate in r:
      i = max_d + mate['d_err']
      j = min(mate['read'].mapping_quality, max_MQ)
      if mate['read_info'].v_list:
        for v_size in mate['read_info'].v_list:
          k = min(max_vlen, max(0, max_vlen + v_size))
          dmv_mat[i, j, k] += 1
      else:
        k = 2 * max_vlen + 1
        dmv_mat[i, j, k] += 1
      dmv_mat[i, j, -1] += 1  # The exact marginal

    yield r


def tlen_matrix():
  pass


@cytoolz.curry
def to_df(riter, tags=None):
  """This is a terminal, it produces a data frame

  :param riter:
  :param tags:
  :return:
  """
  qname = []
  read_data = [
    OrderedDict(
      [
        ('mate', []),
        ('chrom', []),
        ('pos', []),
        ('cigar', []),
        ('MQ', []),
        ('d_err', []),
        ('correct_chrom', []),
        ('correct_pos', []),
        ('correct_cigar', []),
      ] + [
        (tag, [])
        for tag in tags
      ])
    for _ in [0, 1]
  ]

  for r in riter:
    qname.append(r[0]['read'].qname.split('|')[0])
    for n, mate in enumerate(r):
      rdta = read_data[n]
      rd = mate['read']
      rdta['mate'].append(1 if rd.is_read1 else 2)
      rdta['chrom'].append(rd.reference_name)
      rdta['pos'].append(rd.pos)
      rdta['cigar'].append(rd.cigarstring)
      rdta['MQ'].append(rd.mapping_quality)
      rdta['d_err'].append(mate['d_err'])
      rdta['correct_chrom'].append(mate['read_info'].chrom)
      rdta['correct_pos'].append(mate['read_info'].pos)
      rdta['correct_cigar'].append(mate['read_info'].cigar)
      for tag in tags:
        rdta[tag] = rd.get_tag(tag)

  mates_were_paired = True if len(read_data[1]['mate']) else False
  data = OrderedDict(
    [
      (('qname',) if mates_were_paired else 'qname', qname),
    ] + [
      (('mate1', k) if mates_were_paired else k, v)
      for k, v in read_data[0].items()
    ] + ([
      (('mate2', k), v)
      for k, v in read_data[1].items()
    ] if mates_were_paired else [])
  )

  return pd.DataFrame(data)

