import time
from itertools import chain

import pysam
import numpy as np

from mitty.lib.bedfile import read_bed
from mitty.lib.sanitizeseq import sanitize

import logging

logger = logging.getLogger(__name__)


SEED_MAX = (1 << 32) - 1  # Used for seeding rng


def main(fp_out, fasta_fname, sample_name, bed_fname, seed, p_het, models):
  """

  :param fp_out: output file pointer
  :param fasta_fname:
  :param sample_name:
  :param bed_fname:
  :param seed:
  :param p_het:
  :param models:
  :return:
  """
  logger.debug('Starting variant simulation ...')
  t0 = time.time()
  v_cnt_tot = 0

  fp_out.write(generate_header(fasta_fname, sample_name))
  fasta = pysam.FastaFile(fasta_fname)
  rng = np.random.RandomState(seed=seed)
  for region in read_bed(bed_fname):
    t1 = time.time()
    v_cnt = write_out_variants(
      fp_out,
      region,
      [
        model_dispatch[model[0]](
          rng, region=region, seq=fasta.fetch(reference=region[0], start=region[1], end=region[2]),
          p=model[1],
          p_het=p_het,
          min_size=model[2], max_size=model[3])
        for model in models
      ])
    t2 = time.time()
    logger.debug('Wrote {} variants in region {} in {:0.2f}s'.format(v_cnt, region, t2 - t1))
    v_cnt_tot += v_cnt

  t2 = time.time()
  logger.debug('Simulated {} variants in {:0.2f}s'.format(v_cnt_tot, t2 - t0))


def generate_header(fasta_fname, sample_name):
  fasta = pysam.FastaFile(fasta_fname)
  return '\n'.join(
    ['##fileformat=VCFv4.2'] + \
    ['##contig=<ID={}, length={}>'.format(i, l) for i, l in zip(fasta.references, fasta.lengths)] + \
    ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype string">',
     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(sample_name)]) + '\n'


def place_poisson_seq(rng, p, seq):
  """Given a random number generator, a probability and an end point, generate poisson distributed events. Skip bases
  that are 'N'.  For short end_p this may, by chance, generate fewer locations that normal"""
  if p == 0.0:
    return np.array([], dtype='i4')
  l = len(seq)
  return np.array(
    [loc for loc in rng.geometric(p=p, size=int(l * p * 1.2)).cumsum()
     if loc < l and seq[loc] != 'N'], dtype='i4')


def genotype(p_het, rng, l):
  """

  :param rng:
  :param l:
  :return: gt - 0 = 1|0
                1 = 0|1
                2 = 1|1
  """
  r1 = rng.rand(l)  # For the het/hom determination
  r2 = rng.rand(l)  # For hets, for determining which copy
  gt = np.zeros(l, dtype=int)
  gt[r1 > p_het] = 2
  gt[(r1 <= p_het) & (r2 < 0.5)] = 1
  return gt


def snp_model(rng, region, seq, p, p_het, **args):
  """Places SNPs"""
  base_sub = {
    'A': 'CTG',
    'C': 'ATG',
    'T': 'ACG',
    'G': 'ACT'
  }
  pos = place_poisson_seq(rng, p, seq)
  ref = [seq[x] for x in pos]
  alt = [base_sub[r][i] for r, i in zip(ref, rng.randint(0, 3, size=pos.shape[0]))]
  gt = genotype(p_het, rng, pos.shape[0])

  return [pos + region[1] + 1, ref, alt, gt]


t_mat = {
  'A': [0.32654629, 0.17292732, 0.24524503, 0.25528135],
  'T': [0.3489394, 0.25942695, 0.04942584, 0.3422078],
  'G': [0.28778188, 0.21087004, 0.25963262, 0.24171546],
  'C': [0.21644706, 0.20588717, 0.24978216, 0.32788362]
}

cum_t_mat = {
  k: np.cumsum(v)
  for k, v in t_mat.items()
}


def markov_chain(ref, rng, l):
  """

  :param ref:
  :param rng:
  :param l:
  :return:
  """
  dna = 'ACTG'
  alt = [ref] * (l + 1)
  for n, r in enumerate(rng.rand(l)):
    v = cum_t_mat[alt[n]]
    for m in range(4):
      if r <= v[m]:
        alt[n + 1] = dna[m]
        break
    else:
      alt[n + 1] = dna[3]
  return ''.join(alt)


def ins_model(rng, region, seq, p, p_het, min_size, max_size):
  """Insertions uniformly spanning minimum and maximum lengths. Sequences are generated using a Markov chain
generator"""
  pos = place_poisson_seq(rng, p, seq)
  ref = [seq[x] for x in pos]
  alt = [markov_chain(r, rng, l) for r, l in zip(ref, rng.randint(min_size, max_size, size=pos.shape[0]))]
  gt = genotype(p_het, rng, pos.shape[0])
  return [pos + region[1] + 1, ref, alt, gt]


def del_model(rng, region, seq, p, p_het, min_size, max_size):
  """Deletions uniformly spanning minimum and maximum lengths"""
  pos = place_poisson_seq(rng, p, seq)
  ref = [seq[x:x + l + 1] for x, l in zip(pos, rng.randint(min_size, max_size, size=pos.shape[0]))]
  alt = [seq[x] for x in pos]
  gt = genotype(p_het, rng, pos.shape[0])
  return [pos + region[1] + 1, ref, alt, gt]


def copy_ins_model(rng, region, seq, p, p_het, min_size, max_size):
  """The `CINS` model works just like the `INS` model except the insertion sequences, instead of being
novel DNA sequences created with a Markov chain generator, are exact copies of random parts of
the input reference genome. This creates insertions that are more challenging to align to and
assemble, especially when their lengths start to exceed the template size of the sequencing
technology used.
"""
  seq = sanitize(seq)
  pos = place_poisson_seq(rng, p, seq)
  ref = [seq[x] for x in pos]
  alt = [copied_insertion(r, rng, seq, l) for r, l in zip(ref, rng.randint(min_size, max_size, size=pos.shape[0]))]
  gt = genotype(p_het, rng, pos.shape[0])
  try:
    pos, ref, alt, gt = zip(*filter_none_alt(pos + region[1] + 1, ref, alt, gt))
  except ValueError:
    pos, ref, alt, gt = [[], [], [], []]
  return pos, ref, alt, gt


def copied_insertion(ref, rng, seq, l):
  if len(seq) <= l:
    return None
  n0 = rng.randint(len(seq) - l)
  n1 = n0 + l

  if 'N' in seq[n0:n1]:
    return None

  return ref + seq[n0:n1]


def filter_none_alt(pos, ref, alt, gt):
  for p, r, a, g in zip(pos, ref, alt, gt):
    if a is not None:
      yield p, r, a, g


model_dispatch = {
  'SNP': snp_model,
  'INS': ins_model,
  'CINS': copy_ins_model,
  'DEL': del_model
}


def write_out_variants(fp_out, region, v_l):
  """Given a list of list of variants (as returned from the variant model functions) concatenate them and then
   write them out in position order

  :param fp_out:
  :param region:
  :param v_l:
  :return:
  """
  gt_str = ['0|1', '1|0', '1|1']
  contig = region[0]

  pos = list(chain(*(v[0] for v in v_l)))
  ref = list(chain(*(v[1] for v in v_l)))
  alt = list(chain(*(v[2] for v in v_l)))
  gt = list(chain(*(v[3] for v in v_l)))

  for idx in np.argsort(np.array(pos)):  # Sort by position
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT\n
    fp_out.write('{}\t{}\t.\t{}\t{}\t100\tPASS\t.\tGT\t{}\n'.format(contig, pos[idx], ref[idx], alt[idx], gt_str[gt[idx]]))
  return len(pos)