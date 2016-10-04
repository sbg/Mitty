"""Given a FASTQ file generate a table of Base Quality scores"""
from multiprocessing import Process, Queue
import pickle
import time
import logging

import numpy as np
import pysam


logger = logging.getLogger(__name__)
__process_stop_code__ = 'SETECASTRONOMY'


# For debugging
def base_quality_single_threaded(fastq, pkl):
  """Given a fastq file, read through the file

  :param fastq: name of fastq file
  :param pkl: name of pickle file to write to
  :return:
  """
  #                 bp, phred
  bq_mat = np.zeros((500, 100), dtype=np.uint64)
  idx = list(range(500))
  for r in pysam.FastqFile(filename=fastq, persist=False):
    bqa = r.get_quality_array()
    bq_mat[idx[:len(bqa)], bqa] += 1

  pickle.dump(bq_mat, open(pkl, 'wb'))

  return bq_mat


def base_quality(fastq1, fastq2, pkl, threads=2):
  """Given a pair of FASTQs work out the BQ profile for each mate

  :param fastq1: Pair1 and
  :param fastq2: Pair2 fastqs
  :param pkl: name of pickle file to write to
  :param threads: How many 'threads' to use
  :return:
  """
  bq_mats = [
    process_one_fastq(fastq1, threads=threads),
    process_one_fastq(fastq2, threads=threads)
  ]

  pickle.dump(bq_mats, open(pkl, 'wb'))

  return bq_mats


def process_one_fastq(fastq, threads=2, max_reads_in_queue=int(30e6)):
  """

  :param fastq:
  :param threads:
  :param max_reads_in_queue: The default is about 6GB, considering 200 bytes per qual string
  :return:
  """

  t0 = time.time()
  in_queue, out_queue = Queue(max_reads_in_queue), Queue()

  # Start worker processes
  logger.debug('Starting {} threads'.format(threads))
  p_list = [Process(target=process_worker, args=(i, in_queue, out_queue)) for i in range(threads)]
  for p in p_list:
    p.start()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  for r in pysam.FastqFile(filename=fastq, persist=False):
    bqa = r.get_quality_array()
    in_queue.put(bqa)  # TODO: Block when queue gets too big

  # Tell child processes to stop
  logger.debug('Telling child processes to stop')
  for i in range(threads):
    in_queue.put(__process_stop_code__)

  # Get results and add them
  logger.debug('Summing up result matrices')
  bq_mat = out_queue.get()
  for i in range(threads - 1):
    bq_mat += out_queue.get()

  # Wait for workers to finish
  logger.debug('Waiting for workers to shutdown')
  for p in p_list:
    p.join()

  t1 = time.time()
  logger.debug('Finished processing FASTQ in {} s'.format(t1 - t0))

  return bq_mat


def process_worker(worker_no, in_queue, out_queue):
  """Process templates as they are distributed. This is designed to be a worker process for a multiprocessing
  pool, but can be tested without recourse to multiprocessing

  :param worker_no: an id for the worker, not really used in computation
  :param in_queue:  an object with a get method that returns templates when called with next()
  :param out_queue: a queue to put the results on when done
  :return:
  """
  logger.debug('Worker {} starting'.format(worker_no))
  #                 bp, phred
  bq_mat = np.zeros((500, 100), dtype=np.uint64)
  idx = list(range(500))
  for bqa in iter(in_queue.get, __process_stop_code__):
    bq_mat[idx[:len(bqa)], bqa] += 1

  out_queue.put(bq_mat)
  logger.debug('Worker {} stopping'.format(worker_no))