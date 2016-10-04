"""Calculates BQ distribution from a BAM.

The task is parallelized by chromosome.
"""
from multiprocessing import Pool
import time
import pickle
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pysam

logger = logging.getLogger(__name__)


def bq_over_chromosome(bam_fname, chrom):
  """Return BQ matrix for all reads in this chromosome

  :param bam_fname:
  :param chrom: The chromosome id as stored in the BAM header
  :return: 500x100 (read base x qual) uint64 array containing counts
  """
  logger.debug('Processing {}'.format(chrom))
  bq_mat = np.zeros((500, 100), dtype=np.uint64)
  bam_fp = pysam.AlignmentFile(bam_fname, mode='rb')
  for r in bam_fp.fetch(reference=chrom):
    if r.flag > 255: continue  # Ignore secondary/split alignments
    bq_r = r.query_qualities[::-1] if r.flag & 0x10 else r.query_qualities  # Take care of reverse complement
    bq_mat[range(len(bq_r)), bq_r] += 1

  return bq_mat


def process_bam_parallel(bam_fname, pkl, threads=4):
  p = Pool(threads)
  t0 = time.time()
  bam_fp = pysam.AlignmentFile(bam_fname, mode='rb')
  max_chroms = min(24, len(bam_fp.header['SQ']))
  bq_data = {'seq_info': bam_fp.header['SQ']}
  for chrom_data in p.imap_unordered(
    process_bam_section_w,
    ({"bam_fname": bam_fname, "chrom": bam_fp.header['SQ'][i]['SN']}
     for i in range(0, max_chroms))):
    bq_data.update({chrom_data[0]: chrom_data[1]})
    t1 = time.time()
    logger.debug('Processed {} ({} s)'.format(chrom_data[0], t1 - t0))

  pickle.dump(bq_data, open(pkl, 'wb'))

  return bq_data


def process_bam_section_w(args):
  """A thin wrapper to allow proper tracebacks when things go wrong in a thread

  :param args:
  :return:
  """
  import traceback
  try:
    return args['chrom'], bq_over_chromosome(**args)
  except Exception as e:
    traceback.print_exc()
    print('')
    raise e


from matplotlib.colors import LogNorm


def plot_bq_panel(bq_data):
  # TODO: hardcoded for 24 chromosomes, make this data aware?
  seq = bq_data['seq_info']
  fig = plt.figure(figsize=(12, 8))
  rows, cols = 4, 6
  # gc_lim = [0.0, 1.0]
  # cov_lim = [0.0, max_cov]

  for row in range(rows):
    for col in range(cols):
      if row * cols + col > len(seq) - 1: break
      sn = seq[row * cols + col]['SN']
      ax = plt.subplot(rows, cols, row * cols + col + 1)
      plot_bq_metrics(ax, bq_data[sn])
      plt.title(sn)
      # plt.setp(ax, xlim=gc_lim, ylim=cov_lim)

      if row == rows - 1 and col == 0:
        ax.set_xlabel('Read bp')
        ax.set_ylabel('BQ')
      else:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


def plot_bq_combined(bq_data):
  ax = plt.subplot(111)
  plot_bq_metrics(ax, np.sum(bq_mat for k, bq_mat in bq_data.items() if k != 'seq_info'))
  plt.title('Combined')
  ax.set_xlabel('Read bp')
  ax.set_ylabel('BQ')


def plot_bq_metrics(ax, score):
  read_count = score.sum(axis=1)[0]
  max_rlen = score.sum(axis=1).nonzero()[0][-1] + 1
  ax.pcolormesh(score[:max_rlen, :].T, cmap=plt.cm.gray_r)
  ax.plot(np.arange(max_rlen) + 0.5, np.dot(score, np.arange(100))[:max_rlen] / float(read_count), 'y')
  ax.xaxis.set_ticks_position('bottom')
  plt.setp(ax, xlim=(-20, max_rlen + 20), ylim=(0, 60))  # score.shape[1]))


# if __name__ == '__main__':
#   logging.basicConfig(level=logging.DEBUG)
#   bam_fname = '/encrypted/data/notebooks/kghose/bqgc/Projects/bcafc5cc-991c-47db-94f0-da890b680356/NA12878.gathered.bam'
#   process_bam_parallel(bam_fname, 'test.pkl', threads=8)
