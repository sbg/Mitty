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
  d_err = max_d
  if r.is_unmapped:
    d_err = max_d + 2
  elif r.reference_name != ri.chrom:
    d_err = max_d + 1
  else:
    correct_pos = ri.pos
    # Both the strict case and the case where the read comes from inside a long insertion are scored
    # similarly
    if not (strict or ri.cigar[0] == '>'):
      # Score read considering soft clipping
      cigar_op = r.cigartuples[0]
      if cigar_op[0] in [1, 2, 3, 4, 5]:
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
     ('XM', ri.cigar, 'Z'),
     ('XV', ri.v_list, 'B')])