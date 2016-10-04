"""Computes GC content vs coverage from a BAM.

The reference is split up into sections say 10 kb long. For each chunk we compute the GC content and
the average coverage. We return this as an array spanning the reference with GC and coverage values.
This array can be further processed to get metrics useful for modeling the GC bias.

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


def gc_and_coverage_for_region(bam_fp, fasta_fp, region):
  """Compute GC content (as a fraction) and average coverage for this region

  :param bam_fp:
  :param fasta_fp:
  :param region: e.g. "1:10000-20000"
  :return: gc, cov
  """
  seq = fasta_fp.fetch(region=region)
  if seq.count('N') / float(len(seq)) > 0.1:
    return None, None  # Probably centromere or telomere

  gc = float(seq.count('G') + seq.count('C')) / len(seq)
  cov_l = [b.n for b in bam_fp.pileup(region=region)]
  cov = float(sum(cov_l)) / max(len(cov_l), 1)

  return gc, cov


def gc_and_coverage_for_chromosome(bam_fname, fasta_fname, chrom_idx, block_len=10000):
  """

  :param bam_fname:  Passing file names rather than file objects for parallelization
  :param fasta_fname:
  :param chrom_idx: 0, 1, 2 ... referring to chroms in the bam header
  :param block_len: how many bp to chunk by
  :return: an array of gc and cov values
  """
  bam_fp = pysam.AlignmentFile(bam_fname, mode='rb')
  region_start, region_end = 1, bam_fp.header['SQ'][chrom_idx]['LN']
  logger.debug('Processing {}:{}-{}'.format(bam_fp.header['SQ'][chrom_idx]['SN'], region_start, region_end))
  fasta_fp = pysam.FastaFile(fasta_fname)
  return np.array([
    gc_and_coverage_for_region(bam_fp, fasta_fp, region='{}:{}-{}'.format(bam_fp.header['SQ'][chrom_idx]['SN'], r, r + block_len))
    for r in range(region_start, region_end, block_len)
  ], dtype=[('gc', float), ('coverage', float)])


def process_bam_parallel(bam_fname, fasta_fname, pkl, block_len=10000, threads=4):
  p = Pool(threads)
  t0 = time.time()
  bam_fp = pysam.AlignmentFile(bam_fname, mode='rb')
  max_chroms = min(24, len(bam_fp.header['SQ']))
  gc_cov = {'block_len': block_len, 'seq_info': bam_fp.header['SQ']}
  for chrom_data in p.imap_unordered(
    process_bam_section_w,
    ({"bam_fname": bam_fname, "fasta_fname": fasta_fname, "chrom_idx": i, "block_len": block_len}
     for i in range(0, max_chroms))):
    gc_cov.update({bam_fp.header['SQ'][chrom_data[0]]['SN']: chrom_data[1]})
    t1 = time.time()
    logger.debug('Processed {} ({} s)'.format(bam_fp.header['SQ'][chrom_data[0]], t1 - t0))

  pickle.dump(gc_cov, open(pkl, 'wb'))

  return gc_cov


def process_bam_section_w(args):
  """A thin wrapper to allow proper tracebacks when things go wrong in a thread

  :param args:
  :return:
  """
  import traceback
  try:
    # return gc_and_coverage_for_chromosome(bam_fname, fasta_fname, chrom_idx, block_len=10000)
    return (args['chrom_idx'], gc_and_coverage_for_chromosome(**args))
  except Exception as e:
    traceback.print_exc()
    print('')
    raise e


def plot_gc_cov(gc_cov, max_cov=60, title='GC/cov'):
  """Plot a multi panel plot fo GC/cov data

  :param gc_cov:
  :return:
  """
  from matplotlib.colors import LogNorm
  # TODO: hardcoded for 24 chromosomes, make this data aware?
  seq = gc_cov['seq_info']
  fig = plt.figure(figsize=(12, 8))
  rows, cols = 4, 6
  gc_lim = [0.0, 1.0]
  cov_lim = [0.0, max_cov]

  for row in range(rows):
    for col in range(cols):
      if row * cols + col > len(seq) - 1: break
      sn = seq[row * cols + col]['SN']
      ax = plt.subplot(rows, cols, row * cols + col + 1)
      x, y = gc_cov[sn]['gc'], gc_cov[sn]['coverage']
      x, y = x[~(np.isnan(x) | np.isnan(y))], y[~(np.isnan(x) | np.isnan(y))]
      # plt.plot(x, y, 'k.', ms=0.1)
      ax.hist2d(x, y, bins=71, range=[gc_lim, cov_lim], cmap=plt.cm.gray_r, norm=LogNorm())
      plt.title(sn)
      plt.setp(ax, xlim=gc_lim, ylim=cov_lim)

      if row == rows - 1 and col == 0:
        ax.set_xlabel('GC content')
        ax.set_ylabel('Coverage')
      else:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

  fig.suptitle(title)