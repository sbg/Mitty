"""Contains the logic for handling read model corruption invocation"""
from multiprocessing import Process, Queue
import time

import pysam
import numpy as np

from mitty.simulation.sequencing.writefastq import writer, load_qname_sidecar, parse_qname

import logging
logger = logging.getLogger(__name__)


SEED_MAX = (1 << 32) - 1  # Used for seeding rng
__process_stop_code__ = 'SETECASTRONOMY'


def multi_process(read_module, read_model, fastq1_in, fastq1_out, sidecar_in, sidecar_out,
                  fastq2_in=None, fastq2_out=None, processes=2, seed=7):
  """

  :param read_module:
  :param read_model:
  :param fastq1_in:
  :param fastq1_out:
  :param sidecar_in:
  :param sidecar_out:
  :param fastq2_in:
  :param fastq2_out:
  :param processes:
  :param seed:
  :return:
  """

  long_qname_table = load_qname_sidecar(sidecar_in)

  seed_rng = np.random.RandomState(seed)

  logger.debug('Starting {} workers'.format(processes))
  in_queue, out_queue = Queue(10000), Queue(10000)
  p_list = [Process(target=worker,
                    args=(i, read_module, read_model, long_qname_table, in_queue, out_queue, seed_rng.randint(SEED_MAX)))
            for i in range(processes)]
  for p in p_list:
    p.start()

  logger.debug('Starting writer process')
  wr = Process(target=writer, args=(fastq1_out, sidecar_out, fastq2_out, out_queue))
  wr.start()

  t0 = time.time()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  fastq_l = [pysam.FastxFile(fastq1_in)]
  if fastq2_in is not None: fastq_l += [pysam.FastxFile(fastq2_in)]

  cnt = 0
  for cnt, reads in enumerate(zip(*fastq_l)):
    # [(qname, seq, seq) ... ]
    in_queue.put((reads[0].name,) + tuple(r.sequence for r in reads))
    if cnt % 100000 == 0:
      logger.debug('Read {} templates'.format(cnt))

  logger.debug('Stopping child processes')
  for i in range(processes):
    in_queue.put(__process_stop_code__)
  for p in p_list:
    p.join()

  logger.debug('Stopping writer')
  out_queue.put(__process_stop_code__)
  wr.join()

  t1 = time.time()
  logger.debug('Processed {} templates in {:0.2f}s ({:0.2f} t/s)'.format(cnt, t1 - t0, cnt/(t1 - t0)))


def worker(worker_id, read_module, read_model, long_qname_table, in_queue, out_queue, seed):
  """Worker grabs templates from the in_queue, passes them through the read corrupter and then back to the
  out_queue

  :param worker_id: Just an int to help us identify things if needed
  :param read_model:
  :param long_qname_table: a dict mapping qname index to full qname for qnames > 254 chars in length
  :param in_queue:
  :param out_queue:
  :param seed:
  :return:
  """
  corrupt_rng = np.random.RandomState(seed)
  logger.debug('Starting worker {} ...'.format(worker_id))
  cnt, t0 = -1, time.time()
  for cnt, template in enumerate(iter(in_queue.get, __process_stop_code__)):
    # template - (qname, seq, seq)
    ri = parse_qname(template[0], long_qname_table)
    # read_module.corrupt_template gets and returns data in the form
    # (
    #   index,  - same as the original index from the perfect reads
    #   sample_name,
    #   chrom,
    #   copy,
    #   (
    #     (strand, pos, cigar, (v1,v2,...), MD, seq, qual)
    #     ...  [repeated as for as many reads in this template]
    #   )
    # )
    out_queue.put(read_module.corrupt_template(
      read_model,
      (
        ri[0].index,
        ri[0].sample,
        ri[0].chrom,
        ri[0].cpy,
        ((r.strand, r.pos, r.special_cigar or r.cigar, r.v_list, r.md, seq, None) for r, seq in zip(ri, template[1:]))
      ),
      corrupt_rng))
    if cnt % 100000 == 0:
      t1 = time.time()
      logger.debug('Worker {}: Processed {} templates ({} t/s)'.format(worker_id, cnt + 1, (cnt + 1) / (t1 - t0)))
  t1 = time.time()
  logger.debug('... worker {} processed {} templates ({} t/s)'.format(worker_id, cnt + 1, (cnt + 1) / (t1 - t0)))