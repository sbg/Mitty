"""Contains the logic for handling read model corruption invocation"""
from multiprocessing import Process, Queue
import queue
import time

import pysam
import numpy as np


import logging
logger = logging.getLogger(__name__)


SEED_MAX = (1 << 32) - 1  # Used for seeding rng
__process_stop_code__ = 'SETECASTRONOMY'


def multi_process(read_module, read_model, fastq1_in, fastq1_out, fastq2_in=None, fastq2_out=None, processes=2, seed=7):
  """

  :param read_module:
  :param read_model:
  :param fastq1_in:
  :param fastq1_out:
  :param fastq2_in:
  :param fastq2_out:
  :param processes:
  :param seed:
  :return:
  """
  seed_rng = np.random.RandomState(seed)

  logger.debug('Starting {} workers'.format(processes))
  in_queue, out_queue = Queue(), Queue()
  p_list = [Process(target=worker,
                    args=(i, read_module, read_model, in_queue, out_queue, seed_rng.randint(SEED_MAX)))
            for i in range(processes)]
  for p in p_list:
    p.start()

  logger.debug('Starting writer process')
  wr = Process(target=writer, args=(fastq1_out, fastq2_out, out_queue))
  wr.start()

  t0 = time.time()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  fastq_l = [pysam.FastxFile(fastq1_in)]
  if fastq2_in is not None: fastq_l += [pysam.FastxFile(fastq2_in)]

  cnt = 0
  for cnt, reads in enumerate(zip(*fastq_l)):
    # [(seq, seq) ... ]
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


def worker(worker_id, read_module, read_model, in_queue, out_queue, seed):
  """Worker grabs templates from the in_queue, passes them through the read corrupter and then back to the
  out_queue

  :param worker_id: Just an int to help us identify things if needed
  :param read_model:
  :param in_queue:
  :param out_queue:
  :param seed:
  :return:
  """
  corrupt_rng = np.random.RandomState(seed)
  logger.debug('Starting worker {} ...'.format(worker_id))
  cnt, t0 = 0, time.time()
  for cnt, template in enumerate(iter(in_queue.get, __process_stop_code__)):
    out_queue.put(read_module.corrupt_template(read_model, template, corrupt_rng))
    if cnt % 100000 == 0:
      t1 = time.time()
      logger.debug('Worker {}: Processed {} templates ({} t/s)'.format(worker_id, cnt + 1, (cnt + 1) / (t1 - t0)))
  t1 = time.time()
  logger.debug('... worker {} processed {} templates ({} t/s)'.format(worker_id, cnt + 1, (cnt + 1) / (t1 - t0)))


def writer(fastq1_out, fastq2_out=None, out_queue=None):
  """Write templates to file

  :param fastq1_out:
  :param fastq2_out:
  :param out_queue:
  :return:
  """
  logger.debug('Writing out corrupted reads ...')
  t0 = time.time()

  cnt = 0
  fastq_l = [open(fastq1_out, 'w')]
  if fastq2_out is not None: fastq_l += [open(fastq2_out, 'w')]
  for cnt, template in enumerate(iter(out_queue.get, __process_stop_code__)):
    # Data format is [(qname, seq, bq), ...] (Only one read for SE)
    for fp, r in zip(fastq_l, template):
      fp.write('@{}\n{}\n+\n{}\n'.format(r[0], r[1], r[2]))
  for fp in fastq_l:
    fp.close()

  t1 = time.time()
  logger.debug('... finished writing {} templates in {:0.2f}s ({:0.2f} t/s)'.format(cnt, t1 - t0, cnt/(t1 - t0)))