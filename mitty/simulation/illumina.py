"""Example read model

TODO: Figure out the logistics of having diffrent models and loading GC-bias, sequencing error etc.

"""
import numpy as np


SEED_MAX = (1 << 32) - 1  # Used for seeding rng


def read_model_params(gc_bias=None, rlen=150, tlen=500, tlen_std=30, diploid_coverage=30.0):
  """

  :param gc_bias:
  :param rlen:
  :param tlen:
  :param tlen_std:
  :param diploid_coverage: Coverage we expect from a diploid genome
  :return: dict of params,
           passes -> how many times we should call generate reads for each region
  """
  # coverage = read count * read len / genome len
  # coverage = P * read len P = expected number of reads per base
  # P = p * passes    p = probability of read per pass
  # One template for two reads
  # p = P / (2 * passes)
  # p = coverage / (2 * rlen * passes)
  p = 1.0
  passes = 1
  while p > 0.01:
    passes *= 2
    p = 0.5 * diploid_coverage / (2 * rlen * passes)

  return {
    'diploid_coverage': diploid_coverage,
    'rlen': rlen,
    'tlen': tlen,
    'tlen_std': tlen_std,
    'p': p,
    'passes': passes
  }


def generate_reads(model, region, seed=7):
  """

  :param model: dict as returned by read_model_params
  :param region: ('1', 10001, 103863906) Like the BED file
  :param seed:
  :return:
  """
  if not (0 <= seed <= SEED_MAX):
    raise ValueError('Seed value {} is out of range 0 - {}'.format(seed, SEED_MAX))

  seed_rng = np.random.RandomState(seed)
  tloc_rng, tlen_rng, file_order_rng = \
    [np.random.RandomState(s) for s in seed_rng.randint(SEED_MAX, size=3)]

  return _reads_for_template_in_region(
        model['rlen'], file_order_rng,
        *_templates_for_region(
          region, model['p'], model['rlen'], model['tlen'], model['tlen_std'], tloc_rng, tlen_rng))


def _templates_for_region(region, p, rlen, tlen, tlen_std, tloc_rng, tlen_rng):
  # TODO: GC-bias goes here
  est_block_size = int((region[2] - region[1]) * p * 1.2)
  ts = tloc_rng.geometric(p=p, size=est_block_size).cumsum() + region[1]
  tl = (tlen_rng.randn(ts.size) * tlen_std + tlen).astype(int)
  np.clip(tl, rlen, tlen + 5 * tlen_std, out=tl)
  te = ts + tl
  idx = (te < region[2])
  return ts[idx], te[idx]


def _reads_for_template_in_region(rlen, file_order_rng, ts, te):
  """

  :param rlen: Fixed - Illumina reads
  :param file_order_rng: Decides which one of the pair (F or R) comes first in the file (or goes in file1 of the pair)
  :param ts: Start of template
  :param te: End of template
  :return:
  """
  r0l = np.full(ts.size, rlen, dtype=np.uint32)
  r1l = np.full(ts.size, rlen, dtype=np.uint32)

  r0fo = file_order_rng.randint(2, size=ts.size, dtype='i1')
  r1fo = 1 - r0fo

  r0p = ts
  r1p = te - rlen

  return [
    {
      'file_order': r0fo,
      'pos': r0p,
      'len': r0l
    },
    {
      'file_order': r1fo,
      'pos': r1p,
      'len': r1l
    }
  ]


class ReadLocModel(object):
  def __init__(self, gc_bias=None, rlen=150, tlen=500, tlen_std=30, coverage=30.0, passes=10):
    # This will be responsible for loading GC-bias models, template length models etc. from profile files
    # p is the probability of a template starting on any given base
    # passes is how many time we run the algorithm to get the desired coverage
    self.p = 1.0
    passes = max(passes / 2, 1)
    while self.p > 0.1:
      passes *= 2
      self.p = (coverage / tlen) / passes
    self.passes = int(passes)
    self.rlen = rlen
    self.tlen = tlen
    self.tlen_std = tlen_std





  # def generate_templates(self, bed, seed=7):
  #   """
  #
  #   :param bed: [('1', 10001, 103863906),
  #                ('1', 103913907, 206332221),
  #                ('1', 206482222, 249250621) ... ]
  #   :param references: ('1', '2', '3', ... )
  #   :param lengths: (249250621, 243199373, 198022430, ... )
  #   :param seed: for random number generator, for repeatability
  #   :return: [(pt0, pt1, r0s, r0l, r1s, r1l) ..... ]
  #            pt0 -> template start positon
  #            pt1 -> template end position
  #            r0s -> read0 strand (0, 1)
  #            r0l -> read0 length
  #            r1s -> read1 strand (0, 1)
  #            r1l -> read1 length
  #
  #   The bed file is the place to restrict reads to particular regions
  #   """
  #   if not (0 <= seed <= SEED_MAX):
  #     raise ValueError('Seed value {} is out of range 0 - {}'.format(seed, SEED_MAX))
  #
  #   seed_rng = np.random.RandomState(seed)
  #   tloc_rng, tlen_rng, rstrand_rng, file_order_rng = \
  #     [np.random.RandomState(s) for s in seed_rng.randint(SEED_MAX, size=4)]
  #
  #   return [
  #     {
  #       'region': region,
  #       'reads': self._reads_for_template_in_region(
  #           self.rlen, rstrand_rng, file_order_rng,
  #           *self._templates_for_region(
  #             region, self.p, self.passes, self.rlen, self.tlen, self.tlen_std, tloc_rng, tlen_rng))
  #     }
  #     for region in bed
  #   ]

  def generate_templates(self, bed, seed=7):
    """

    :param bed: [('1', 10001, 103863906),
                 ('1', 103913907, 206332221),
                 ('1', 206482222, 249250621) ... ]
    :param references: ('1', '2', '3', ... )
    :param lengths: (249250621, 243199373, 198022430, ... )
    :param seed: for random number generator, for repeatability
    :return: [(pt0, pt1, r0s, r0l, r1s, r1l) ..... ]
             pt0 -> template start positon
             pt1 -> template end position
             r0s -> read0 strand (0, 1)
             r0l -> read0 length
             r1s -> read1 strand (0, 1)
             r1l -> read1 length

    The bed file is the place to restrict reads to particular regions
    """
    if not (0 <= seed <= SEED_MAX):
      raise ValueError('Seed value {} is out of range 0 - {}'.format(seed, SEED_MAX))

    seed_rng = np.random.RandomState(seed)
    tloc_rng, tlen_rng, rstrand_rng, file_order_rng = \
      [np.random.RandomState(s) for s in seed_rng.randint(SEED_MAX, size=4)]

    return [
      {
        'region': region,
        'reads': self._reads_for_template_in_region(
            self.rlen, rstrand_rng, file_order_rng,
            *self._templates_for_region(
              region, self.p, self.passes, self.rlen, self.tlen, self.tlen_std, tloc_rng, tlen_rng))
      }
      for region in bed
    ]


  @staticmethod
  def _templates_for_region(region, p, passes, rlen, tlen, tlen_std, tloc_rng, tlen_rng):
    # TODO: GC-bias goes here
    est_block_size = int((region[2] - region[1]) * p * 1.2)
    ts = np.concatenate([
      tloc_rng.geometric(p=p, size=est_block_size).cumsum() + region[1]
      for _ in range(passes)])
    ts.sort(kind='mergesort')
    tl = tlen_rng.randn(ts.size) * tlen_std + tlen
    np.clip(tl, rlen, tlen + 5 * tlen_std, out=tl)
    te = ts + tl
    idx = (te < region[2])
    return ts[idx], te[idx]

  @staticmethod
  def _reads_for_template_in_region(rlen, rstrand_rng, file_order_rng, ts, te):
    """

    :param rlen: Fixed - Illumina reads
    :param rstrand_rng: Decides whether read comes from forward or reverse strand
    :param file_order_rng: Decides whether this read goes in first or second place in file
                           Basically means read
    :param ts: Start of template
    :param te: End of template
    :return:
    """
    r0s = rstrand_rng.randint(2, size=ts.size, dtype='i1')
    r1s = 1 - r0s

    r0l = np.full(ts.size, rlen, dtype=np.uint32)
    r1l = np.full(ts.size, rlen, dtype=np.uint32)

    r0fo = file_order_rng.randint(2, size=ts.size, dtype='i1')
    r1fo = 1 - r0fo

    r0p = ts
    r1p = te - rlen

    return [
      {
        'strand': r0s,
        'file_order': r0fo,
        'pos': r0p,
        'len': r0l
      },
      {
        'strand': r1s,
        'file_order': r1fo,
        'pos': r1p,
        'len': r1l
      }
    ]