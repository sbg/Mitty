from multiprocessing import Process, Queue
import time
import logging
import os

import pysam

from mitty.lib.bamfetch import fetch
from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname

__process_stop_code__ = 'SETECASTRONOMY'
logger = logging.getLogger(__name__)


def main(bam_fname, sidecar_fname, out_fname,
         d_range=(-200, 200), reject_d_range=False,
         v_range=(-200, 200), reject_v_range=False,
         reject_reads_with_variants=False,
         reject_reference_reads=False,
         strict_scoring=False, do_not_index=True, processes=2):
  """This function extracts poorly aligned reads from a BAM from simulated reads"""
  # Set up the I/O queues and place all BAM contigs on the work queue
  work_q = Queue()
  for ref in pysam.AlignmentFile(bam_fname).references:
    work_q.put(ref)
  for _ in range(processes):
    work_q.put(__process_stop_code__)

  # Start workers
  long_qname_table = load_qname_sidecar(sidecar_fname)
  file_fragments = ['{}.{:04}.bam'.format(out_fname, i) for i in range(processes)]
  p_list = [
    Process(target=worker,
            args=(i, bam_fname, long_qname_table, file_fragments[i],
                  d_range, reject_d_range,
                  v_range, reject_v_range,
                  reject_reads_with_variants,
                  reject_reference_reads,
                  strict_scoring, work_q))
    for i in range(processes)
  ]
  for p in p_list:
    p.start()

  # Orderly exit
  for p in p_list:
    p.join()

  merge_sorted_fragments(out_fname, file_fragments, do_not_index)


def worker(worker_id, bam_fname, long_qname_table, bam_out,
           d_range, reject_d_range,
           v_range, reject_v_range,
           reject_reads_with_variants,
           reject_reference_reads,
           strict_scoring=False, in_q=None):
  logger.debug('Starting worker {} ...'.format(worker_id))

  t0, tot_cnt = time.time(), 0

  bam_fp = pysam.AlignmentFile(bam_fname)
  out_fp = pysam.AlignmentFile(bam_out, 'wb', header=bam_fp.header)

  for reference in iter(in_q.get, __process_stop_code__):
    logger.debug('Worker {}: Contig {} ...'.format(worker_id, reference))
    t1 = time.time()
    cnt = process_contig(bam_fp, reference, long_qname_table,
                         d_range, reject_d_range,
                         v_range, reject_v_range,
                         reject_reads_with_variants,
                         reject_reference_reads,
                         strict_scoring, out_fp)
    t2 = time.time()
    logger.debug(
      'Worker {}: Contig {}: {} reads in {:2f}s ({:2f} r/s)'.format(worker_id, reference, cnt, t2 - t1, cnt / (t2 - t1)))
    tot_cnt += cnt

  out_fp.close()
  t1 = time.time()
  logger.debug(
    'Worker {}: Processed {} reads in {:2f}s ({:2f} r/s)'.format(worker_id, tot_cnt, t1 - t0, tot_cnt / (t1 - t0)))

  logger.debug('Sorting {} -> {}'.format(bam_out, bam_out + '.sorted'))
  t0 = time.time()
  pysam.sort('-m', '1G', '-o', bam_out + '.sorted', bam_out)
  os.remove(bam_out)
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))


def process_contig(bam_fp, reference, long_qname_table,
                   d_range, reject_d_range,
                   v_range, reject_v_range,
                   reject_reads_with_variants,
                   reject_reference_reads,
                   strict_scoring, out_fp):
  cnt = 0
  max_d = d_range[1] + 10000

  for r in fetch(bam_fp, reference=reference):
    if r.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    cnt += 1
    ri = parse_qname(r.qname, long_qname_table=long_qname_table)[1 if r.is_read2 else 0]
    d_err = score_alignment_error(r, ri=ri, max_d=max_d, strict=strict_scoring)

    is_ref_read = len(ri.v_list) == 0
    if is_ref_read and reject_reference_reads:
      continue

    if not is_ref_read and reject_reads_with_variants:
      continue

    in_d_err_range = d_range[0] <= d_err <= d_range[1]
    if in_d_err_range == reject_d_range:
      continue

    if not is_ref_read:
      keep_read = False
      for v in ri.v_list:
        in_v_range = v_range[0] < v < v_range[1]
        if in_v_range != reject_v_range:
          keep_read = True
          break
      if not keep_read:
        continue

    r.set_tag('XD', d_err)
    out_fp.write(r)

  return cnt


def merge_sorted_fragments(bam_fname, file_fragments, do_not_index):
  logger.debug('Merging sorted BAM fragments ...')
  t0 = time.time()
  pysam.merge('-rpcf', bam_fname, *[f + '.sorted' for f in file_fragments])
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('Removing fragments')
  for f in file_fragments:
   os.remove(f + '.sorted')

  if not do_not_index:
    logger.debug('BAM index {} ...'.format(bam_fname))
    t0 = time.time()
    pysam.index(bam_fname, bam_fname + '.bai')
    t1 = time.time()
    logger.debug('... {:0.2f}s'.format(t1 - t0))