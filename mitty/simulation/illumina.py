"""Example read model

TODO: Figure out the logistics of having diffrent models and loading GC-bias, sequencing error etc.

"""
import numpy as np


SEED_MAX = (1 << 32) - 1  # Used for seeding rng


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
    chrom_cpy_rng, tloc_rng, tlen_rng, rstrand_rng = \
      [np.random.RandomState(s) for s in seed_rng.randint(SEED_MAX, size=4)]

    return [
      {
        'region': region,
        'reads': self._reads_for_template_in_region(
            self.rlen, chrom_cpy_rng, rstrand_rng,
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
  def _reads_for_template_in_region(rlen, chrom_cpy_rng, rstrand_rng, ts, te):
    cc = chrom_cpy_rng.randint(1, 3, size=ts.size, dtype='i1')
    r0s = rstrand_rng.randint(2, size=ts.size, dtype='i1')
    r1s = 1 - r0s
    r0l = np.full(ts.size, rlen, dtype=np.uint32)
    r1l = np.full(ts.size, rlen, dtype=np.uint32)

    return np.rec.fromarrays(
      (ts, te, cc, r0s, r0l, r1s, r1l),
      dtype=[
        ('tplt_st', 'u4'),
        ('tplt_end', 'u4'),
        ('chrom_cpy', 'i1'),
        ('r0_strand', 'i1'),
        ('r0_len', 'u4'),
        ('r1_strand', 'i1'),
        ('r1_len', 'u4')])