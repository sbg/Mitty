"""
Algorithm (made interesting by having to deal with paired reads)

- Consider paired reads (or single reads) as a unit (a list).
- If any one of the unit passes the read filter, write out the whole unit
- For unpaired reads, this is as far as the algorithm needs to go

- For paired reads in the same chromosome, this is easily done by keeping the read/read-pair in a dictionary
  until both pairs are found and then flushing it out
- For reads whose mates are in different chromosomes or one mate is unmapped we have some fun
  - Two workers will run into this qname (since it appears in two chromosomes)
  - Our data can't cross the process boundary
- After we have finished processing a chromosome
  - we send the qname, contig and pos of all residual reads back to the parent


OK - let's try a single threaded version since everything I can think of with multiple
processes eventually requires us to go through the BAM again for the residual reads since
reads can't cross process boundaries.

The single threaded version would simply keep everything in a dict until flushed out.
As a bonus, the code is considerably simpler for the single process version
"""
import time
import logging
import os

import pysam

from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname

__process_stop_code__ = 'SETECASTRONOMY'
logger = logging.getLogger(__name__)


def main(bam_fname, sidecar_fname, out_fname,
         d_range=(-200, 200), reject_d_range=False,
         v_range=(-200, 200), reject_v_range=False,
         reject_reads_with_variants=False,
         reject_reference_reads=False,
         strict_scoring=False, do_not_index=True, processes=2):
  """This function extracts reads from a simulation BAM that match the filter critera

  :param bam_fname:
  :param sidecar_fname:
  :param out_fname:
  :param d_range:
  :param reject_d_range:
  :param v_range:
  :param reject_v_range:
  :param reject_reads_with_variants:
  :param reject_reference_reads:
  :param strict_scoring:
  :param do_not_index:
  :param processes:
  :return:
  """
  def _filter_pass(_r):
    """

    :param _r:
    :return: T/F, d_err
    """
    ri = parse_qname(_r.qname, long_qname_table=long_qname_table)[1 if _r.is_read2 else 0]

    is_ref_read = len(ri.v_list) == 0
    if is_ref_read and reject_reference_reads:
      return False, 0

    if not is_ref_read and reject_reads_with_variants:
      return False, 0

    _d_err = score_alignment_error(_r, ri=ri, max_d=max_d, strict=strict_scoring)

    in_d_err_range = d_range[0] <= _d_err <= d_range[1]
    if in_d_err_range == reject_d_range:
      return False, 0

    if not is_ref_read:
      # All variants are inside/outside v_range and we want to/do not want to reject the range
      if all((v_range[0] <= v <= v_range[1]) == reject_v_range for v in ri.v_list):
        return False, 0

    return True, _d_err

  se_bam = is_single_end_bam(bam_fname)
  bam_fp = pysam.AlignmentFile(bam_fname)
  long_qname_table = load_qname_sidecar(sidecar_fname)

  unsorted_out_fname = out_fname + '.unsorted'
  out_fp = pysam.AlignmentFile(unsorted_out_fname, 'wb', header=bam_fp.header)

  in_cnt = 0
  max_d = d_range[1] + 10000
  read_dict = {}

  t0 = time.time()
  for rd in bam_fp.fetch(until_eof=True):
    if rd.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    in_cnt += 1
    if in_cnt % 1000000 == 0:
      t1 = time.time()
      logger.debug(
        'Processed {} reads in {:2f}s ({:2f} r/s) {}'.format(
          in_cnt, t1 - t0, in_cnt / (t1 - t0), '' if se_bam else '(dict size {})'.format(len(read_dict))))

    if se_bam:
      keep, d_err = _filter_pass(rd)
      if keep:
        rd.set_tag('XD', d_err)
        out_fp.write(rd)
    else:
      if rd.qname[:20] not in read_dict:
        read_dict[rd.qname[:20]] = [None, None]

      rl = read_dict[rd.qname[:20]]
      rl[0 if rd.is_read1 else 1] = rd

      if all(rl):
        keep1, d_err1 = _filter_pass(rl[0])
        keep2, d_err2 = _filter_pass(rl[1])
        if keep1 or keep2:
          rl[0].set_tag('XD', d_err1)
          rl[1].set_tag('XD', d_err2)
          out_fp.write(rl[0])
          out_fp.write(rl[1])
        del read_dict[rd.qname[:20]]

  out_fp.close()
  t1 = time.time()
  logger.debug(
    'Processed {} reads in {:2f}s ({:2f} r/s) {}'.format(
      in_cnt, t1 - t0, in_cnt / (t1 - t0), '' if se_bam else '(dict size {})'.format(len(read_dict))))

  logger.debug('Sorting {} -> {}'.format(unsorted_out_fname, out_fname))
  t0 = time.time()
  pysam.sort('-m', '1G', '-o', out_fname, unsorted_out_fname)
  os.remove(unsorted_out_fname)
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  if not do_not_index:
    logger.debug('BAM index {} ...'.format(bam_fname))
    t0 = time.time()
    pysam.index(out_fname, out_fname + '.bai')
    t1 = time.time()
    logger.debug('... {:0.2f}s'.format(t1 - t0))


def is_single_end_bam(bam_fname):
  bam_fp = pysam.AlignmentFile(bam_fname)
  r = next(bam_fp, None)
  return not r.is_paired if r is not None else True  # Empty BAM? Don't care