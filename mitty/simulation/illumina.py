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


def generate_reads(model, p_min, p_max, seed=7):
  """

  :param model: dict as returned by read_model_params
  :param p_min: Start of reads
  :param p_max: End of reads
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
          p_min, p_max, model['p'], model['rlen'], model['tlen'], model['tlen_std'], tloc_rng, tlen_rng))


def _templates_for_region(p_min, p_max, p, rlen, tlen, tlen_std, tloc_rng, tlen_rng):
  # TODO: GC-bias goes here
  est_block_size = int((p_max - p_min) * p * 1.2)
  ts = tloc_rng.geometric(p=p, size=est_block_size).cumsum() + p_min + 1
  tl = (tlen_rng.randn(ts.shape[0]) * tlen_std + tlen).astype(int)
  np.clip(tl, rlen, tlen + 5 * tlen_std, out=tl)
  te = ts + tl
  idx = (te < p_max)
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