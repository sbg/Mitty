"""Code for assessing alignment accuracy"""
import re

import pysam

from mitty.simulation.sequencing.writefastq import ri, load_qname_sidecar, parse_qname
# good to have these together. Code that uses score_alignment_error is very likely to use
# ri, load_qname_sidecar and parse_qname

cigar_parser = re.compile(r'(\d+)(\D)')


def score_alignment_error(r, ri, max_d=200, strict=False):
  """
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
    if strict or ri.cigar[0] == '>':
      # Both the strict case and the case where the read comes from inside a long insertion are scored
      # similarly
      d_err = max(min((r.pos + 1 - ri.pos), max_d), -max_d)
    else:  # Score read as correct if read is placed at any breakpoint correctly
      correct_pos = ri.pos
      for cnt, op in cigar_parser.findall(ri.cigar):
        if op == '=' or op == 'M' or op == 'X':
          if abs(r.pos + 1 - correct_pos) < abs(d_err):
            d_err = r.pos + 1 - correct_pos
          correct_pos += int(cnt)
        elif op == 'D':
          correct_pos += int(cnt)

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


def main(bam_in, bam_out, max_d=200, strict=False, sample_name=None, processes=2, sidecar_fname=None, truth_bam=None):
  """If simulated reads, sidecar_fname must be passed (and there can't be any truth BAM). The score is
  written into the Xd tag. If a truth_bam is passed then this bam is used as the reference and Xd is computed with
  respect to this.

  :param bam_in:
  :param bam_out:
  :param max_d:
  :param strict:
  :param sample_name:
  :param processes:
  :param sidecar_fname:
  :param truth_bam:
  :return:
  """
  pass


def d_err(bam_in, bam_out, max_d=200, strict=False, sidecar_fname=None):
  """

  :param bam_in:
  :param bam_out:
  :param max_d:
  :param strict:
  :param sidecar_fname:
  :return:
  """
  fp = pysam.AlignmentFile(bam_in, mode='rb')
  long_qname_table = load_qname_sidecar(sidecar_fname) if sidecar_fname else None
  fp_out = pysam.AlignmentFile(bam_out, mode='wb', header=fp.header)
  for r in fp.fetch(until_eof=True):  # This includes the unmapped reads
    ri = parse_qname(r.qname, long_qname_table)[0 if r.is_read1 else 1]
    tag_alignment(r, ri, max_d, strict)
    fp_out.write(r)


def p_diff(bam_in, bam_out, max_d=200, strict=False, truth_bam=None):
  fp = pysam.AlignmentFile(bam_in, mode='rb')
  fp_truth = pysam.AlignmentFile(truth_bam, mode='rb')





  fp_out = pysam.AlignmentFile(bam_out, mode='wb', header=fp.header)
  for r in fp.fetch(until_eof=True):  # This includes the unmapped reads
    ri = parse_qname(r.qname, long_qname_table)[0 if r.is_read1 else 1]
    tag_alignment(r, ri, max_d, strict)
    fp_out.write(r)

def iterate_over_bams(bam_fp_l):
  pass
