"""This function rips through a BAM from simulated reads and bins reads into a three dimensional histogram.

The dimensions are:
 Xd - alignment error
 MQ - mapping quality
 vlen - length of variant carried by read

The histogram is saved as a numpy array and different views of the histogram are plotted.
"""
from multiprocessing import Pool
import time
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pysam
import numpy as np

#from mitty.simulation.sequencing.writefastq import load_qname_sidecar, parse_qname
from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname

logger = logging.getLogger(__name__)


def main(bam_fname, sidecar_fname, max_xd=200, max_MQ=70, strict_scoring=False,  max_vlen=200, processes=2):
  """This function rips through a BAM from simulated reads and bins reads into a three dimensional histogram.

  The dimensions are:
    Xd - alignment error  [0]  -max_xd, ... 0, ... +max_xd, wrong_chrom, unmapped (2 * max_xd + 3)
    MQ - mapping quality  [1]  0, ... max_MQ (max_MQ + 1)
    vlen - length of variant carried by read [2]  Ref, < -max_vlen , -max_vlen, ... 0, ... +max_vlen, > +max_vlen
                                                  ( 2 * max_vlen + 1 + 2 + 1)

  :param bam_fname:
  :param sidecar_fname:
  :param max_xd:
  :param max_MQ:
  :param strict_scoring:
  :param max_vlen:
  :param processes:
  :return:
  """
  long_qname_table = load_qname_sidecar(sidecar_fname)
  xmv_mat = None
  pool = Pool(processes=processes)
  for this_xmv_mat in pool.imap_unordered(
    worker, ((bam_fname, long_qname_table, ref, max_xd, max_MQ, max_vlen)
             for ref in pysam.AlignmentFile(bam_fname).references)):
    xmv_mat = this_xmv_mat if xmv_mat is None else xmv_mat + this_xmv_mat
  return xmv_mat


def save(xmv_mat, fname):
  np.save(fname, xmv_mat)


def worker(args):
  return process_fragment(*args)


def process_fragment(bam_fname, long_qname_table, reference, max_xd=200, max_MQ=70, max_vlen=200, strict_scoring=False):
  """

  :param bam_fname:
  :param long_qname_table:
  :param contig:
  :param start:
  :param end:
  :param max_xd:
  :param max_MQ:
  :param max_vlen:
  :param strict_scoring:
  :return:
  """
  logger.debug('Processing {} ...'.format(reference))

  t0, cnt = time.time(), -1

  xmv_mat = np.zeros(shape=(2 * max_xd + 3, max_MQ + 1, 2 * max_vlen + 1 + 2 + 1), dtype=int)
  v_off_idx, max_v_idx = max_vlen + 2, 2 * max_vlen + 3
  bam_fp = pysam.AlignmentFile(bam_fname)
  for cnt, r in enumerate(bam_fp.fetch(reference=reference)):
    ri = parse_qname(r.qname, long_qname_table=long_qname_table)[1 if r.is_read2 else 0]
    d_err = score_alignment_error(r, ri=ri, max_d=max_xd, strict=strict_scoring)
    if not ri.v_list:
      xmv_mat[max_xd + d_err, min(r.mapping_quality, max_MQ), 0] += 1
    else:
      for v_size in ri.v_list:
        xmv_mat[max_xd + d_err, min(r.mapping_quality, max_MQ), min(max_v_idx, max(1, v_off_idx + v_size))] += 1

  t1 = time.time()
  logger.debug('Processed {}: {} reads in {:2f}s ({:2f} r/s)'.format(reference, cnt + 1, t1 - t0, (cnt + 1) / (t1 - t0)))

  return xmv_mat


def plot_panel(xmv_mat, fig_fname=None):
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99)

  ax1 = plt.subplot(411)
  plot_mean_MQ_vs_derr(ax1, xmv_mat, show_ax_label=True)

  if fig_fname is not None:
    plt.savefig(fig_fname)


# Different plots
def plot_mean_MQ_vs_derr(ax, xmv_mat, plot_bin_size=5, show_ax_label=False):
  """Plot mean_MQ (y-axis) against d_err (x-axis) for given range of variant sizes

  :param xmv_mat:
  :param vlen_slice:
  :return:
  """
  mq_mat = xmv_mat.sum(axis=2)
  data_cnt = mq_mat.sum(axis=1)
  mq_vector = np.arange(mq_mat.shape[1])
  mean_mq = np.dot(mq_vector, mq_mat.T) / np.clip(data_cnt, 1, data_cnt.max() + 1).astype(float)

  max_derr = int((mean_mq.shape[0] - 3) / 2)
  x1, x2 = max_derr * 1.5, max_derr * 1.7
  xt = [-max_derr, -int(max_derr/2), 0, int(max_derr/2), max_derr]

  ax.plot(range(-max_derr, max_derr + 1), mean_mq[:2 * max_derr + 1], 'k.')
  ax.plot([x1, x2], mean_mq[2 * max_derr + 1:], 'ko')
  ax.set_xticks(xt + [x1, x2])
  ax.set_xticklabels(xt + ['WC', 'UM'] if show_ax_label else [])


def plot_mq(mq_mat, plot):
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99)

  ax1 = plt.subplot(411)
  plot_mq_vs_p_err(mq_mat, ax1, yscale='linear')

  ax2 = plt.subplot(412, sharex=ax1)
  plot_mq_vs_p_err(mq_mat, ax2, yscale='log')

  ax3 = plt.subplot(413, sharex=ax1)
  plot_mq_vs_read_cnt(mq_mat, ax3)

  ax4 = plt.subplot(414)
  plot_derr_vs_mq(mq_mat, ax4)

  plt.savefig(plot)


def plot_mq_vs_p_err(mq_mat, ax, yscale='log'):
  mq = np.arange(100)
  for fmt, d_margin, label in zip(['g*', 'ro'], [0, 10], [r'$d_{err} = 0$', r'$d_{err} < 10$']):
    p_err = compute_p_err(mq_mat, d_margin)
    ax.plot(mq, p_err, fmt, mfc='none', label=label)

  ax.plot(mq, 10 ** (-mq / 10), 'k:')
  plt.setp(ax, xlabel='MQ', ylabel=r'$p_{err}$', xlim=(-1, 65), yscale=yscale)
  ax.legend(loc='upper right' if yscale=='linear' else 'lower left')


def compute_p_err(mq_mat, d_margin=0):
  d_max = int((mq_mat.shape[0] - 3)/2)
  c = mq_mat[d_max - d_margin:d_max + 1 + d_margin, :, :].sum(axis=0).sum(axis=1)
  t = mq_mat.sum(axis=0).sum(axis=1)
  return (t - c) / np.clip(t, 1, 1 + t.max())


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