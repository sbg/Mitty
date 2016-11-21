import time

import pysam
import numpy as np

import logging
logger = logging.getLogger(__name__)


# NOTE TO SELF: Though it seems like we could convert the psyam variant handling structures to something
# that looks like it is more useful for us re: read generation, it is not worth the effort and the loss
# in modularity and maintainability. It is better to simply have functions that convert the VCF records
# to the structure used by read generation at a later stage.


# The pysam variant class is a little too convoluted for our use, so we simplify the fields to this
# and add some useful information
# Note the absence of the CHROM field
class Variant(object):
  __slots__ = ('pos', 'ref', 'alt', 'cigarop', 'oplen')

  def __init__(self, pos, ref, alt, cigarop, oplen):
    self.pos = pos
    self.ref = ref
    self.alt = alt
    self.cigarop = cigarop
    self.oplen = oplen

  def tuple(self):
    return self.pos, self.ref, self.alt, self.cigarop, self.oplen

  def __repr__(self):
    return self.tuple().__repr__()


# Follows standard BED convention: https://genome.ucsc.edu/FAQ/FAQformat#format1
# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# chromStart - The starting position of the feature in the chromosome or scaffold.
#              The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold.
#            The chromEnd base is not included in the display of the feature.
#
# For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100,
# and span the bases numbered 0-99.
def read_bed(bed_fname):
  return list(map(lambda x: (x[0], int(x[1]), int(x[2])), map(lambda x: x.split(), open(bed_fname, 'r').readlines())))


# Unphased variants always go into chrom copy 0|1
# We get results in the form of a ploid_bed
def load_variant_file(fname, sample, bed_fname):
  """Use pysam to read in a VCF file and convert it into a form suitable for use in Mitty

  :param fname:
  :param sample:
  :return: dict of numpy recarrays
  """
  mode = 'rb' if fname.endswith('bcf') else 'r'
  vcf_fp = pysam.VariantFile(fname, mode)
  vcf_fp.subset_samples([sample])
  return [
    split_copies(region, [v for v in vcf_fp.fetch(contig=region[0], start=region[1], stop=region[2])])
    for region in read_bed(bed_fname)
  ]


def split_copies(region, vl):
  """Given a list of pysam.cbcf.VariantRecord split it into multiple lists, one for each chromosome copy

  :param vl:
  :return: [c1|0, c0|1]
  """
  # Sniff out the ploidy
  if len(vl) == 0:
    logger.warning('Empty region ({}), assuming diploid'.format(region))
    ploidy = 2
  else:
    ploidy = len(vl[0].samples[0]['GT'])
    logger.debug('Region: {}, ploidy: {}'.format(region, ploidy))

  # cpy_l = [
  #   (cpy, '|'.join(['0'] * cpy + ['1'] + ['0'] * (ploidy - 1 - cpy)))
  #   for cpy in range(ploidy)
  # ]
  #
  # return {
  #   'region': region,
  #   'v': dict(
  #     [
  #       (gt, list(filter(None, (parse(v, cpy=cpy) for v in vl))))
  #       for cpy, gt in cpy_l
  #     ]
  #   )
  # }

  return {
    'region': region,
    'v': [
        list(filter(None, (parse(v, cpy=cpy) for v in vl)))
        for cpy in range(ploidy)
    ]
  }


def parse(v, cpy):
  """Take a pysam.cbcf.VariantRecord and convert it into a Variant object

  :param v: variant
  :param cpy: 0 = 1|0, 1 = 0|1 OR
              0 = 1|0|0, 1 = 0|1|0, 2 = 0|0|1 etc
  :return: Variant(object)
  """
  if v.samples[0]['GT'][cpy] == 0:  # Not present in this copy
    return None
  alt = v.samples[0].alleles[cpy]
  l_r, l_a = len(v.ref), len(alt)
  if l_r == 1:
    if l_a == 1:
      op, op_len = 'X', 0
    else:
      op, op_len = 'I', l_a - l_r
  elif l_a == 1:
    op, op_len = 'D', l_r - l_a
  else:
    raise ValueError("Complex variants present in VCF. Please filter or refactor these.")

  return Variant(v.pos, v.ref, v.samples[0].alleles[cpy], op, op_len)


def prepare_variant_file(fname_in, sample, bed_fname, fname_out, write_mode='w'):
  """Prepare a variant file with only the given sample, complex variant calls filtered out,
  and restricted to the given bed file

  :param fname_in:
  :param sample:
  :param bed_fname:
  :param fname_out:
  :return: - output is to file
  """
  def _complex_variant(_v):
    var = _v.samples.values()[0]
    for alt in var.alleles:
      if _v.rlen > 1 and len(alt) > 1:
        if _v.ref != alt:
          logger.debug('Filtered out {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, var.alleles))
          return True
    return False

  logger.debug('Starting filtering ...')
  t0 = time.time()

  mode = 'rb' if fname_in.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(fname_in, mode)
  vcf_in.subset_samples([sample])
  vcf_out = pysam.VariantFile(fname_out, mode=write_mode, header=vcf_in.header)
  v_cnt, fltr_cnt = 0, 0
  for region in read_bed(bed_fname):
    logger.debug('Filtering {}'.format(region))
    n = 0
    for n, v in enumerate(vcf_in.fetch(contig=region[0], start=region[1], stop=region[2])):
      if _complex_variant(v):
        fltr_cnt += 1
        continue
      vcf_out.write(v)
    v_cnt += n

  logger.debug('Processed {} variants'.format(v_cnt))
  logger.debug('Filtered out {} complex variants'.format(fltr_cnt))
  t1 = time.time()
  logger.debug('Took {} s'.format(t1 - t0))