"""
This contains code to create, fill and operate on the 8D alignment histogram and slices thereof
"""
from collections import OrderedDict, Counter
import logging

import cytoolz
import numpy as np
import pysam
import xarray as xr

from mitty.benchmarking.alignmentscore import score_alignment_error, correct_tlen, load_qname_sidecar, parse_qname


logger = logging.getLogger(__name__)


def fastmultihist(sample, bins):
  """N-dimensional histogram

  :param sample:  N x D array (samples x dimensions)
  :param bins:    D length list of arrays representing bin edges
                  Some of D may be np.array(None), in which case those dimensions are ignored
                  and a D' dimensional histogram is returned with the None dimensions
                  collapsed
  :return:

  Because np.histogramdd is broken in terms of performance.
  Because we don't need no nanny state
  """
  N, D = sample.shape

  # Compute the bin number each sample falls into.
  Ncount = [np.digitize(sample[:, i], bins[i]) - 1 for i in np.arange(D) if bins[i].shape != ()]

  # Unravel the D-dimensional coordinates onto a 1-D index
  nbin = np.array([len(b) - 1 for b in bins if b.shape != ()], int)
  xy = np.zeros(N, int)
  for i in np.arange(0, len(nbin) - 1):
    xy += Ncount[i] * nbin[i + 1:].prod()
  xy += Ncount[-1]

  # Histogram the indices and then ravel the index counts back into D-dimensions
  return np.bincount(xy, minlength=nbin.prod()).reshape(nbin)


def make_tick_labels(bin_edges, var_label='x'):
  """Utility function that allows us to set bin labels for our histogram

  :param bin_edges: the bin edges
  :param var_label: what do we want to label this 'x', 'y', 'MQ' etc.
  :return:
  """
  if bin_edges is None: return None
  b = bin_edges
  return [
    '{} <= {} < {}'.format(round(b[i], 2), var_label, round(b[i+1], 2))
    for i in range(len(b) - 1)
  ]


def default_histogram_parameters():
  """Give us reasonable defaults for the full (8D) histogram.

  :return: A dictionary which we can then modify if we wish
           This dictionary can then be passed to `initialize_histogram`
  """
  return {
    'xd_bin_edges': (-np.inf, -100, -40, -20, -10, -5, 0, 1, 6, 11, 21, 41, 101, np.inf),
    'max_d': 200,
    'mq_bin_edges': (0, 1, 11, 21, 31, 41, 51, 61),
    'v_bin_edges': (-np.inf, -20, 0, 1, 21, np.inf),
    'xt_bin_edges': (-np.inf, -100, -40, -20, -10, -5, 0, 1, 6, 11, 21, 41, 101, np.inf),
    't_bin_edges': (-np.inf, 200, 301, 401, 601, np.inf),
    'name': 'PairedAlignmentHistogram'
  }


# Some notes on histogram sizes saved as gzipped pickle file
# (15, 15, 5, 5, 7, 7, 13, 5) -> 150 kB
def initialize_histogram(
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

  **We may set some of the dimensions to None to reduce the dimensionality**
  """
  # Some of the dimensions need extra bins for particular exceptions to the data

  if xd_bin_edges is not None:
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
  else:
    xd_labels = None

  if v_bin_edges is not None:
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
  else:
    max_v = None
    v_labels = None

  # The full suite of dimensions
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

  # This is bin edges in the form of a list so that histogram can use it
  # We make this list before collapsing None dimensions.
  # histograminator needs the None dimensions to know which dimensions of the input data to ignore
  full_bin_edges = [k[2] for k in dimensions]

  # Reduce dimensions, if applicable
  dimensions = [d for d in dimensions if d[3] is not None]

  metadata = {
    'max_d': max_d,
    'max_v': max_v,
    'description': [k[1] for k in dimensions],
    'bin_centers': [(k[2][:-1] + k[2][1:]) / 2 for k in dimensions],
    'bin_edges': full_bin_edges  # Note, not reduced dimensionality, but with 'None' where dimensions are collapsed
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
def histograminator(titer, histogram_def=None, buf_size=1000000):
  """Iterate over read pairs and fill out the histogram

  :param titer:
  :param histogram_def: This can be a list of histogram defs. If there are more than one
                        Multiple histograms will be filled out. The constraint is that the
                        max_d and max_v of each histogram have to match
  :param buf_size:
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
    """given a buffer of data bin it into the histogram
    """
    for pah in aah:
      pah.data += fastmultihist(buf[:idx, :], pah.attrs['bin_edges'])

  assert histogram_def is not None, 'Histogram definition must be passed'
  if not isinstance(histogram_def, (list, tuple)):
    histogram_def = (histogram_def, )

  aah = [initialize_histogram(**hd) for hd in histogram_def]
  max_d_set = set(a.attrs['max_d'] for a in aah)
  max_v_set = set(a.attrs['max_v'] for a in aah if a.attrs['max_v'] is not None)
  max_d = max_d_set.pop()
  max_v = max_v_set.pop() or 1  # We could have a case where we never bin by the max_v

  if len(max_d_set) > 0:
    raise RuntimeError('max_d for each histogram has to be the same')

  if len(max_v_set) > 0:
    raise RuntimeError('max_v for each histogram has to be the same')

  temp_read = pysam.AlignedSegment()

  bs = buf_size
  buf = np.empty(shape=(bs, 8))
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

  return aah


def save_histogram(pah, fname):
  """Saves gzipped to pickle. Nice compression, imperceptibly slower than unzipped version but
  a LOT smaller. xarray's native serialization formats are hillariously broken.

  :param pah:
  :param fname:
  :return:
  """
  # from a hint http://henrysmac.org/blog/2010/3/15/python-pickle-example-including-gzip-for-compression.html
  with gzip.open(fname, 'wb') as f:
    pickle.dump(pah, f, protocol=-1)


def load_histogram(fname):
  with gzip.open(fname, 'rb') as f:
    return pickle.load(f)


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

  new_pah = pah[slices].sum(axis=collapse_idx_l)
  new_pah.attrs = pah.attrs  # sum does not copy these over, unfortunately
  return new_pah


def force_sum(p1, p2):
  """p1 + p2 when the axes labels are different results in an inner product (which is strange)
  this function forces an element-wise addition and then applies axes metadata from p1 to the result

  :param p1:
  :param p2:
  :return:
  """
  return xr.DataArray(
    data=p1.data + p2.data,
    coords=p1.coords,
    dims=p1.dims,
    name=p1.name,
    attrs=p1.attrs
  )
