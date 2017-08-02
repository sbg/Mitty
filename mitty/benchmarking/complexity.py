"""
The function vcf_complexity goes through a VCF annotating variants with local complexity measures. It is
slow - the C++ translation of the program is much faster.
"""

import pickle
import time
import logging
from bisect import bisect_left

import pysam
import numpy as np

logger = logging.getLogger(__name__)

DNA = 'ACGT'


# Code for generating test cases
def de_bruijn(k, n):
  """
  de Bruijn sequence for alphabet k
  and subsequences of length n.
  """
  alphabet = k
  k = len(k)

  a = [0] * k * n
  sequence = []

  def db(t, p):
    if t > n:
      if n % p == 0:
        sequence.extend(a[1:p + 1])
    else:
      a[t] = a[t - p]
      db(t + 1, p)
      for j in range(a[t - p] + 1, k):
        a[t] = j
        db(t + 1, t)

  db(1, 1)
  sequence.extend(sequence[:n - 1])

  return "".join(alphabet[i] for i in sequence)


def dna_alphabet(k=1):
  word = [0] * (k + 1)
  while word[-1] is 0:
    yield ''.join(DNA[w] for w in word[:-1])
    for n in range(k + 1):
      word[n] += 1
      if word[n] < 4:
        break
      else:
        word[n] = 0


class ShannonEntropy:
  """Computes Shannon Entropy of k-mers in given sequence. Initialize with value of k.
  Shannon entropy of a sequence ranges from >0 to 4**k)
  """
  def __init__(self, k=1):
    self.k = k
    self.alphabet = list(dna_alphabet(k=k))

  def header_info(self):
    """Return a list of elements that should be added to the VCF header so:

    vcf_in.header.info.add(*se.header_info())
    """
    return 'SE{}'.format(self.k), 1, 'Float', '{}-mer Shannon entropy for the local region around the variant'.format(self.k)

  def complexity(self, s, **kwargs):
    l = len(s) - self.k + 1
    p = [s.count(a) / l for a in self.alphabet]
    return sum(- px * np.log2(max(px, 1e-6)) for px in p)


class LinguisticComplexity:
  """Computes linguistic complexity"""
  @staticmethod
  def header_info():
    """Return a list of elements that should be added to the VCF header so:

    vcf_in.header.info.add(*lc.header_info())
    """
    return 'LC', 1, 'Float', 'Linguistic complexity of local sequence'

  @staticmethod
  def complexity(s, **kwargs):
    """Compute the linguistic complexity of sequence s

    :param s:
    :return:
    """
    num, den = 1, 1
    for k in range(1, len(s)):
      k4 = 4**k # For DNA
      num += min(len(set(s[i:i+k] for i in range(len(s) - k + 1))), k4)
      den += min(len(s) - k + 1, k4)
    return num / den


# Code for converting a BedGraph file into a numpy data - structure
# ./bigWigToBedGraph wgEncodeCrgMapabilityAlign100mer.bigWig map100.bg && gzip map100.bg
#
def read_bedgraph(bg):
  current_chrom, data = None, []
  with open(bg, 'r') as fp:
    for ln in fp:
      chrom, _, stop, val = ln.split('\t')
      if chrom != current_chrom:
        if current_chrom is not None:
          yield current_chrom, np.array(data, dtype=[('pos', np.uint32), ('mappability', np.float16)])
        current_chrom, data = chrom, []
      data += [(int(stop), float(val[:-1]))]

  if current_chrom is not None:
    yield current_chrom, np.array(data, dtype=[('pos', np.uint32), ('mappability', np.float16)])


# map100.bg (~682MB) takes 62.026s to load
# map100.bg.gz is ~126 MB
# map100.bg.pkl (~130 MB) takes ~ 1s to load
def bedgraph2numpy(bg):
  return {
    chrom: data
    for chrom, data in read_bedgraph(bg)
  }


class Mappability:
  """Attaches a mappability score to a variant"""
  def __init__(self, bed_graph=None):
    data = pickle.load(open(bed_graph, 'rb'))
    self.pos_data = {k: data[k]['pos'] for k in data.keys()}
    self.map_data = {k: data[k]['mappability'] for k in data.keys()}

  @staticmethod
  def header_info():
    """Return a list of elements that should be added to the VCF header so:

    vcf_in.header.info.add(*mp.header_info())
    """
    return 'MAP', 1, 'Float', 'Mappability at variant location'

  def complexity(self, s, chrom=None, pos=None):
    # return float(self.data[chrom][np.searchsorted(self.data[chrom]['pos'], pos, side='left')]['mappability'])
    # https://stackoverflow.com/questions/15139299/performance-of-numpy-searchsorted-is-poor-on-structured-arrays
    return float(self.map_data[chrom][bisect_left(self.pos_data[chrom], pos)])


# This is painfully slow ...
# LC alone runs at 400 v/s
# SE1 alone runs at 1382 v/s
# SE2 alone runs at 1310 v/s
# SE4 alone runs at 834.87 v/s (54226.4 v/s with C++ code. 65x Yeah baby!)
# MAP alone runs at 120 v/s with np.searchsorted and at 1400 v/s with bisect
# Everything together 170 v/s
def vcf_complexity(vcf_in_fname, vcf_out_fname, ref_fname, bg_fname=None, window_size=100):
  """

  :param vcf_in_fname:
  :param vcf_out_fname:
  :param ref_fname:
  :param bg_fname:
  :param d:
  :return:
  """
  measures = [ShannonEntropy(k=k) for k in [1, 2, 4]] + [LinguisticComplexity(), Mappability(bed_graph=bg_fname)]
  info_tags = [m.header_info()[0] for m in measures]

  logger.debug('Starting annotation ...')
  t0 = time.time()

  mode = 'rb' if vcf_in_fname.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(vcf_in_fname, mode)
  for m in measures:
    vcf_in.header.info.add(*m.header_info())
  vcf_out = pysam.VariantFile(vcf_out_fname, mode='w', header=vcf_in.header)

  ref_fp = pysam.FastaFile(ref_fname)

  d_2 = int(window_size / 2)
  v_cnt = 0
  for v_cnt, v in enumerate(vcf_in):
    for tag, m in zip(info_tags, measures):
      v.info[tag] = m.complexity(
        s=ref_fp.fetch(reference=v.contig, start=v.start - d_2, end=v.stop + d_2).upper(),
        chrom=v.contig, pos=v.pos)
    vcf_out.write(v)
    if v_cnt % 10000 == 9999:
      t1 = time.time()
      logger.debug('Processed {} variants in {:0.2f}s ({:0.2f}v/s)'.format(v_cnt + 1, t1 - t0, (v_cnt + 1)/(t1 - t0)))

  t1 = time.time()
  logger.debug('Processed {} variants in {:0.2f}s ({:0.2f}v/s)'.format(v_cnt + 1, t1 - t0, (v_cnt + 1) / (t1 - t0)))
