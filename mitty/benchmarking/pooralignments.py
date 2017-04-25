from multiprocessing import Process, Queue
import time
import logging
import os

import pysam

from mitty.lib.bamfetch import fetch
from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname

__process_stop_code__ = 'SETECASTRONOMY'
logger = logging.getLogger(__name__)


def main(bam_fname, sidecar_fname, out_fname, xd_threshold=1, strict_scoring=False, sort_and_index=True, processes=2):
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
    Process(target=worker, args=(i, bam_fname, long_qname_table, file_fragments[i], xd_threshold, strict_scoring, work_q, sort_and_index))
    for i in range(processes)
  ]
  for p in p_list:
    p.start()

  # Orderly exit
  for p in p_list:
    p.join()

  if sort_and_index:
    merge_sorted_fragments(out_fname, file_fragments)


def worker(worker_id, bam_fname, long_qname_table, bam_out, xd_threshold, strict_scoring=False, in_q=None, sort=True):
  """

  :param worker_id:
  :param bam_fname:
  :param long_qname_table:
  :param bam_out:
  :param xd_threshold:
  :param strict_scoring:
  :param in_q:
  :return:
  """
  logger.debug('Starting worker {} ...'.format(worker_id))

  t0, tot_cnt = time.time(), 0

  bam_fp = pysam.AlignmentFile(bam_fname)
  out_fp = pysam.AlignmentFile(bam_out, 'wb', header=bam_fp.header)

  for reference in iter(in_q.get, __process_stop_code__):
    logger.debug('Worker {}: Contig {} ...'.format(worker_id, reference))
    t1 = time.time()
    cnt = process_contig(bam_fp, reference, long_qname_table, xd_threshold, strict_scoring, out_fp)
    t2 = time.time()
    logger.debug(
      'Worker {}: Contig {}: {} reads in {:2f}s ({:2f} r/s)'.format(worker_id, reference, cnt, t2 - t1, cnt / (t2 - t1)))
    tot_cnt += cnt

  out_fp.close()
  t1 = time.time()
  logger.debug(
    'Worker {}: Processed {} reads in {:2f}s ({:2f} r/s)'.format(worker_id, tot_cnt, t1 - t0, tot_cnt / (t1 - t0)))

  if sort:
    logger.debug('Sorting {} -> {}'.format(bam_out, bam_out + '.sorted'))
    t0 = time.time()
    pysam.sort('-m', '1G', '-o', bam_out + '.sorted', bam_out)
    os.remove(bam_out)
    t1 = time.time()
    logger.debug('... {:0.2f}s'.format(t1 - t0))


def process_contig(bam_fp, reference, long_qname_table, xd_threshold, strict_scoring, out_fp):
  cnt = 0
  max_d = xd_threshold + 10000
  for r in fetch(bam_fp, reference=reference):
    if r.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    cnt += 1
    ri = parse_qname(r.qname, long_qname_table=long_qname_table)[1 if r.is_read2 else 0]
    d_err = score_alignment_error(r, ri=ri, max_d=max_d, strict=strict_scoring)
    if abs(d_err) >= xd_threshold:
      r.set_tag('XD', d_err)
      out_fp.write(r)

  return cnt


def merge_sorted_fragments(bam_fname, file_fragments):
  logger.debug('Merging sorted BAM fragments ...')
  t0 = time.time()
  pysam.merge('-rpcf', bam_fname, *[f + '.sorted' for f in file_fragments])
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('Removing fragments')
  for f in file_fragments:
   os.remove(f + '.sorted')

  logger.debug('BAM index {} ...'.format(bam_fname))
  t0 = time.time()
  pysam.index(bam_fname)
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))