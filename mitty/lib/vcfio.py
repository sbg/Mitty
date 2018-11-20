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
  v_check = UnusableVariantFilter(ploidy)
  return list(filter(None, (parse(v, cpy=cpy, v_check=v_check) for v in vl)))


def parse(v, cpy, v_check):
  """Take a pysam.cbcf.VariantRecord and convert it into a Variant object

  :param v: variant
  :param cpy: 0 = 1|0, 1 = 0|1 OR
              0 = 1|0|0, 1 = 0|1|0, 2 = 0|0|1 etc
  :return: Variant(object)
  """
  if v.samples[0]['GT'][cpy] == 0:  # Not present in this copy
    return None

  if v_check.unusable(v):
    logger.error("Unusable variants present in VCF. Please filter or refactor these.")
    raise ValueError("Unusable variants present in VCF. Please filter or refactor these.")

  alt = v.samples[0].alleles[cpy]
  if alt:
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
  else:
      return None


class UnusableVariantFilter:
  def __init__(self, ploidy):
    self.p_overlap = [0] * ploidy
    self.last_variant = [(0, '', '') for _ in range(ploidy)]

  def unusable(self, v):
    """Return True if we can use this variant

    :param v:
    :return:
    """
    var = v.samples.values()[0]
    is_unusable = \
      self._complex_variant(v) or \
      self._angle_bracketed_id(v, var) or \
      self._breakend_replacement(v) or \
      self._illegal_overlap(v, var, self.p_overlap) or \
      self._duplicate_variant(v, var, self.last_variant)

    if not is_unusable:
      for n, (g, alt) in enumerate(zip(var['GT'], var.alleles)):
        # lv -> (pos, ref, alt)
        if g:
          self.p_overlap[n] = v.stop - 1
          self.last_variant[n] = (v.pos, v.ref, alt)

    return is_unusable

  @staticmethod
  def _complex_variant(_v):
    for alt in _v.alts:
      if _v.rlen > 1 and len(alt) > 1:
        logger.debug('Complex variant {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, _v.alts))
        return True
    return False

  @staticmethod
  def _angle_bracketed_id(_v, var):
    for alt in var.alleles:
      if alt and alt[0] == '<':
        logger.debug('Angle bracketed variant entry {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, var.alleles))
        return True
    return False

  @staticmethod
  def _breakend_replacement(_v):
    if _v.info.get('SVTYPE', None) == 'BND':
      logger.debug('Breakend entry {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, _v.alts))
      return True
    return False

  @staticmethod
  def _illegal_overlap(_v, var, _p_overlap):
    is_illegal = False
    for n, (g, alt, po) in enumerate(zip(var['GT'], var.alleles, _p_overlap)):
      # _v.start and po are 0-indexed
      if g:
        if len(alt) == len(_v.ref) == 1:  # SNP
          start = _v.start  # A SNP's footprint is where it is
        else:
          start = _v.start + 1 # INS and DEL affect one base over
        if start <= po:  # This is overlapping don't use the variant
          is_illegal = True
          logger.debug('Illegal overlap {}:{} {} -> {} (previous variant ends at {})'.format(_v.contig, _v.pos, _v.ref, _v.alts, po + 1))
          break

    return is_illegal

  @staticmethod
  def _duplicate_variant(_v, var, _last_variant):
    is_duplicate = False
    for n, (g, alt, lv) in enumerate(zip(var['GT'], var.alleles, _last_variant)):
      # lv -> (pos, ref, alt)
      if g:
        if (lv[0] == _v.pos) & (lv[1] == _v.ref) & (lv[2] == alt):
          is_duplicate = True
          logger.debug(
            'Duplicate line {}:{} {} -> {}'.format(_v.contig, _v.pos, _v.ref, _v.alts))
          break
    return is_duplicate


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

    v_check = UnusableVariantFilter(sniff_ploidy(vcf_in, contig=region[0]))

    for n, v in enumerate(vcf_in.fetch(contig=region[0], start=region[1], stop=region[2])):
      if not any(v.samples.values()[0]['GT']): continue  # This variant does not exist in this sample
      if v_check.unusable(v):
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