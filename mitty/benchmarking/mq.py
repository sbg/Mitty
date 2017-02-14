"""Go through BAM file made from simulated reads. Score alignment of each read and create a data structure
consisting of MQ and simple d metric. Plot this data"""
from multiprocessing import Pool
import time
import logging

import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm

from mitty.benchmarking.alignmentscore import load_qname_sidecar, parse_qname, score_alignment_error


logger = logging.getLogger(__name__)


def process_bam(bam_fname, max_d=200, strict=False, sample_name=None, processes=2, out_csv=None):
  bam_fp = pysam.AlignmentFile(bam_fname)

  p = Pool(processes=processes)
  mq_mat = np.zeros((100, 2 * max_d + 1), dtype=np.uint64)
  t_st = t0 = time.time()
  tot_rd = 0
  for shard in p.imap_unordered(process_bam_part_w,
                                ((bam_fname, s['SN'], max_d, strict, sample_name) for s in bam_fp.header['SQ'])):
    rd_cnt = shard.sum()
    tot_rd += rd_cnt
    mq_mat += shard
    dt = time.time() - t0
    logger.debug('{} templates processed ({:0.2f} templates/sec).'.format(rd_cnt, rd_cnt / dt))
    t0 = time.time()

  t_nd = time.time()
  logger.debug('{} templates processed in {:0.2f}s ({:0.2f} templates/sec).'.format(tot_rd, t_nd - t_st, tot_rd / (t_nd - t_st)))

  if out_csv is not None:
    logger.debug('Saving MQ matrix to {}'.format(out_csv))
    np.savetxt(out_csv, mq_mat, delimiter=',', fmt='%d')

  return mq_mat


def process_bam_part_w(args):
  return process_bam_part(*args)


def process_bam_part(bam_fname, reference, max_d=200, strict=False, sample_name=None):
  """Work out MQ and D for

  :param bam_fname:
  :param reference:
  :param max_d:
  :param strict:
  :param sample_name:
  :return:
  """
  mq_mat = np.zeros((100, 2 * max_d + 1), dtype=np.uint64)
  bam_fp = pysam.AlignmentFile(bam_fname)
  for r in bam_fp.fetch(reference=reference):
    if sample_name is not None and not r.qname.startswith(sample_name): continue
    ri = parse_qname(r.qname)[1 if r.is_read2 else 0]  # TODO: Will this work for SE reads?
    d_err = score_alignment_error(r, ri, max_d=max_d,
                                  strict=strict)  # d_err = max(min((r.pos + 1 - ri.pos if r.reference_name == ri.chrom else max_d), max_d), -max_d)
    mq_mat[r.mapq, d_err + max_d] += 1
  return mq_mat


def plot_mq(mq_mat, plot):
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99)

  ax1 = plt.subplot(411)
  plot_mq_vs_perr(mq_mat, ax1, yscale='linear')

  ax2 = plt.subplot(412, sharex=ax1)
  plot_mq_vs_perr(mq_mat, ax2, yscale='log')

  ax3 = plt.subplot(413, sharex=ax1)
  plot_mq_vs_read_cnt(mq_mat, ax3)

  ax4 = plt.subplot(414)
  plot_derr_vs_mq(mq_mat, ax4)

  plt.savefig(plot)


def plot_mq_vs_perr(mq_mat, ax, yscale='log'):
  mq = np.arange(100)
  for fmt, d_margin, label in zip(['g*', 'ro'], [0, 10], [r'$d_{err} = 0$', r'$d_{err} < 10$']):
    p_err = compute_p_err(mq_mat, d_margin)
    ax.plot(mq, p_err, fmt, mfc='none', label=label)

  ax.plot(mq, 10 ** (-mq / 10), 'k:')
  plt.setp(ax, xlabel='MQ', ylabel=r'$p_{err}$', xlim=(-1, 65), yscale=yscale)
  ax.legend(loc='upper right' if yscale=='linear' else 'lower left')


def compute_p_err(mq_mat, d_margin=0):
  d_max = int((mq_mat.shape[1] - 1) / 2)
  c = mq_mat[:, d_max - d_margin:d_max + 1 + d_margin].sum(axis=1)
  w = mq_mat[:, :d_max - d_margin].sum(axis=1) + mq_mat[:, d_max + 1 + d_margin:].sum(axis=1)
  return w / np.clip((c + w), 1, (c + w).max())


def plot_derr_vs_mq(mq_mat, ax):
  d_max = (mq_mat.shape[1] - 1) / 2
  d = np.arange(-d_max, d_max + 1)
  cnt = mq_mat.sum(axis=0)
  mq = np.arange(mq_mat.shape[0])
  mean_mq = np.dot(mq, mq_mat) / np.clip(cnt, 1, cnt.max()).astype(float)
  ax.plot(d, mean_mq, 'k.')
  ax.vlines(0, 0, 100, linestyles=':')
  plt.setp(ax, xlabel=r'$d_{err}$', ylabel='MQ (mean)', xlim=(-d_max - 1, d_max + 1))


def plot_mq_vs_read_cnt(mq_mat, ax):
  mq = np.arange(100)
  cnt = mq_mat.sum(axis=1)
  ax.step(mq, cnt, 'k', where='mid')
  plt.setp(ax, xlabel='MQ', ylabel=r'Read count', xlim=(-1, 65), yscale='log')

