"""Code for assessing alignment accuracy"""
import re
import pysam
from mitty.simulation.sequencing.writefastq import ri, load_qname_sidecar, parse_qname
# good to have these together. Code that uses score_alignment_error is very likely to use
# ri, load_qname_sidecar and parse_qname
from mitty.lib.cigars import cigarv2_v1


cigar_parser = re.compile(r'(\d+)(\D)')


def is_simple_case(cigar1, cigar2):
  _, _, op1, _ = cigar_parser.split(cigar1, maxsplit=1)
  _, _, op2, _ = cigar_parser.split(cigar2, maxsplit=1)
  return op1 in ['=', 'M', 'X'] and op2 in ['M', 'X', '=']


def find_first_common_reference_matching_base_positions(r1, r2):
  return next(
    filter(
      lambda x: x[0] is not None and x[1] is not None,
      zip(r1.get_reference_positions(full_length=True),
          r2.get_reference_positions(full_length=True))),
    (None, None))


def score_alignment_error(r, ri, max_d=200, strict=False):
  """Algorithm:

  If strict is True: Look at the distance between the simulated P and the aligned P
  If strict is False:
   If from inside an insertion, treat like strict
   Else, find the first base in the read that is placed on the reference for both the
         aligned and correct reads and compute the difference
         If there is no such base, treat like strict

  Look at the aligned CIGAR and compute the difference between exact P and
                     aligned P after accounting for any soft clip at the start of the read

  :param r: aligned read
  :param ri: readinfo for correct alignment
  :param strict: If True, simply compute difference between simulated P and aligned P
                 if False, find first common reference matching base
  :return: -max_d <= d_err <= max_d + 2
           d_err = max_d + 1 if wrong chrom
           d_err = max_d + 2 if unmapped
  """
  if r.is_unmapped:
    d_err = max_d + 2
  elif r.reference_name != ri.chrom:
    d_err = max_d + 1
  else:
    # This is what we score against when we are 'strict' or inside an insertion
    # Or we can't find a common reference matching base in the correct and aligned reads
    correct_pos, aligned_pos = ri.pos - 1, r.pos
    # If we are not strict AND not in an insertion we use first common reference matching base
    if not strict and ri.special_cigar is None and not is_simple_case(ri.cigar, r.cigarstring):
      rc = pysam.AlignedSegment()
      rc.pos, rc.cigarstring = ri.pos - 1, cigarv2_v1(ri.cigar)
      _correct_pos, _aligned_pos = find_first_common_reference_matching_base_positions(rc, r)
      if _correct_pos is not None:
        correct_pos, aligned_pos = _correct_pos, _aligned_pos

    d_err = max(min((aligned_pos - correct_pos), max_d), -max_d)

  return d_err


def tag_alignment(r, ri, max_d=200, strict=False):
  """Given correct alignment set tags on the read indicating correct alignment and other metadata

  :param r:
  :param ri:
  :return: mutates r
  """
  r.set_tags(
    [('Xd', score_alignment_error(r, ri, max_d, strict), 'i'),
     ('XR', ri.chrom, 'Z'),
     ('XP', ri.pos, 'i'),
     ('XM', ri.cigar, 'Z')] + ([('XV', ri.v_list)] if ri.v_list else []))