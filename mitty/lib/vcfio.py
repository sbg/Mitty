import logging
import time

import pysam

from mitty.lib.bedfile import read_bed

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
    split_copies(region,
                 [v for v in vcf_fp.fetch(contig=region[0], start=region[1], stop=region[2])],
                 sniff_ploidy(vcf_fp, region[0]))
    for region in read_bed(bed_fname)
    ]


def sniff_ploidy(vcf_fp, contig):
  v = next(vcf_fp.fetch(contig=contig), None)
  ploidy = 2
  if v is not None:
    ploidy = len(v.samples[0]['GT'])
  logger.debug(
    'Contig: {}, ploidy: {} {}'.format(contig, ploidy,
                                       '(Assumed. Contig was empty)' if v is None else ''))
  return ploidy


def fetch_first_variant_in_contig_as_empty(vcf_fp, contig):
  v = next(vcf_fp.fetch(contig=contig), None)
  if v is not None:
    v.samples[0]['GT'] = (0,) * len(v.samples[0]['GT'])
  return v


def split_copies(region, vl, ploidy):
  """Given a list of pysam.cbcf.VariantRecord split it into multiple lists, one for each chromosome copy

  :param vl:
  :return: [c1|0, c0|1]
  """
  return {
    'region': region,
    'v': [
        parse_vl(vl, cpy=cpy, ploidy=ploidy)
        for cpy in range(ploidy)
    ]
  }


def parse_vl(vl, cpy, ploidy):
  unusable_variant.p_overlap = [0] * ploidy
  return list(filter(None, (parse(v, cpy=cpy) for v in vl)))


def parse(v, cpy):
  """Take a pysam.cbcf.VariantRecord and convert it into a Variant object

  :param v: variant
  :param cpy: 0 = 1|0, 1 = 0|1 OR
              0 = 1|0|0, 1 = 0|1|0, 2 = 0|0|1 etc
  :return: Variant(object)
  """
  if v.samples[0]['GT'][cpy] == 0:  # Not present in this copy
    return None

  if unusable_variant(v):
    logger.error("Unusable variants present in VCF. Please filter or refactor these.")
    raise ValueError("Unusable variants present in VCF. Please filter or refactor these.")

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


def unusable_variant(v):
  """

  :param v:
  :param p_overlap: meant to be a private mutable (updated in place)
                    Quintisomy should be enough for everyone
  :return: T/F

  This function uses a state variable p_overlap which should be initialized to [0, 0, 0, ....]
  when calling this on a new region or other stretch of variants. The number of elements you
  keep must be at least the ploidy of the genomes you are handling.
  """
  def _complex_variant(_v):
    for alt in _v.alts:
      if _v.rlen > 1 and len(alt) > 1:
        logger.debug('Complex variant {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, _v.alts))
        return True
    return False

  def _angle_bracketed_id(_v):
    var = _v.samples.values()[0]
    for alt in var.alleles:
      if alt[0] == '<':
        logger.debug('Angle bracketed variant entry {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, var.alleles))
        return True
    return False

  def _breakend_replacement(_v):
    if _v.info.get('SVTYPE', None) == 'BND':
      logger.debug('Breakend entry {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, _v.alts))
      return True
    return False

  def _illegal_overlap(_v, _p_overlap):
    is_illegal = False
    var = _v.samples.values()[0]
    for g, po in zip(var['GT'], _p_overlap):
      if g and po > _v.start:  # stop is 1 past the position, for some reason
        is_illegal = True
        logger.debug('Illegal overlap {}:{} {} -> {} (previous variant ends at {})'.format(_v.contig, _v.pos, _v.ref, _v.alts, po))
        break
    else:  # Only gets here if is_illegal is false
      for n, g in enumerate(var['GT']):
        if g:
          _p_overlap[n] = _v.stop
    return is_illegal

  return _complex_variant(v) or _angle_bracketed_id(v) or _breakend_replacement(v) or _illegal_overlap(v, unusable_variant.p_overlap)


def prepare_variant_file(fname_in, sample, bed_fname, fname_out, write_mode='w'):
  """Prepare a variant file with only the given sample, complex and illegal variant calls
  filtered out, and restricted to the given bed file

  :param fname_in:
  :param sample:
  :param bed_fname:
  :param fname_out:
  :return: - output is to file
  """
  logger.debug('Starting filtering ...')
  t0 = time.time()

  mode = 'rb' if fname_in.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(fname_in, mode)
  vcf_in.subset_samples([sample])
  vcf_out = pysam.VariantFile(fname_out, mode=write_mode, header=vcf_in.header)
  processed_cnt, exclude_cnt, include_cnt = 0, 0, 0
  contig_dict = set()

  for region in read_bed(bed_fname):
    logger.debug('Filtering {}'.format(region))
    n, this_include_cnt = -1, 0
    empty_gt = None
    if region[0] not in contig_dict:
      empty_gt = fetch_first_variant_in_contig_as_empty(vcf_in, region[0])
      contig_dict.add(region[0])
    unusable_variant.p_overlap = [0] * sniff_ploidy(vcf_in, contig=region[0])

    for n, v in enumerate(vcf_in.fetch(contig=region[0], start=region[1], stop=region[2])):
      if not any(v.samples.values()[0]['GT']): continue  # This variant does not exist in this sample
      if unusable_variant(v):
        exclude_cnt += 1
        continue
      vcf_out.write(v)
      this_include_cnt += 1

    if this_include_cnt == 0 and empty_gt is not None:
      vcf_out.write(empty_gt)

    include_cnt += this_include_cnt
    processed_cnt += (n + 1)

  logger.debug('Processed {} variants'.format(processed_cnt))
  logger.debug('Sample had {} variants'.format(exclude_cnt + include_cnt))
  logger.debug('Discarded {} variants'.format(exclude_cnt))
  t1 = time.time()
  logger.debug('Took {} s'.format(t1 - t0))