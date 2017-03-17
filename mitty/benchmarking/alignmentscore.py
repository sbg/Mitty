"""Code for assessing alignment accuracy"""
import re

import pysam

from mitty.simulation.sequencing.writefastq import ri, load_qname_sidecar, parse_qname
# good to have these together. Code that uses score_alignment_error is very likely to use
# ri, load_qname_sidecar and parse_qname

cigar_parser = re.compile(r'(\d+)(\D)')


def score_alignment_error(r, ri, max_d=200, strict=False):
  """Algorithm:

  If strict is True: Look at the distance between the simulated P and the aligned P
  If strict is False: Look at the aligned CIGAR and compute the difference between exact P and
                     aligned P after accounting for any soft clip at the start of the read

  :param r: aligned read
  :param ri: readinfo for correct alignment
  :param strict: If True, soft clipped alignments or split alignments are marked incorrect
                 if False alignment to breakpoints gets full score
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
    correct_pos = ri.pos
    # If we are not strict AND not in an insertion we account for left side soft-clipping
    if not strict and ri.special_cigar is None:
      # Score read considering soft clipping
      cigar_op = r.cigartuples[0]
      if 0 < cigar_op[0] < 6:
        correct_pos += cigar_op[1]

    d_err = max(min((r.pos + 1 - correct_pos), max_d), -max_d)

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