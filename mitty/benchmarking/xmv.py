"""This function rips through a BAM from simulated reads and bins reads into a three dimensional histogram.

The dimensions are:
 Xd - alignment error
 MQ - mapping quality
 vlen - length of variant carried by read

The histogram is saved as a numpy array and different views of the histogram are plotted.
"""
from multiprocessing import Process, Queue
import time
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pysam
import numpy as np

from mitty.lib.bamfetch import fetch
from mitty.benchmarking.alignmentscore import score_alignment_error, load_qname_sidecar, parse_qname
from mitty.benchmarking.plot.byvsize import plot_panels

__process_stop_code__ = 'SETECASTRONOMY'
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
  # Set up the I/O queues and place all BAM contigs on the work queue
  work_q, result_q = Queue(), Queue()
  for ref in pysam.AlignmentFile(bam_fname).references:
    work_q.put(ref)
  for _ in range(processes):
    work_q.put(__process_stop_code__)

  # Start workers
  long_qname_table = load_qname_sidecar(sidecar_fname)
  p_list = [
    Process(target=worker, args=(i, bam_fname, long_qname_table, max_xd, max_MQ, max_vlen, strict_scoring,
                                 work_q, result_q))
    for i in range(processes)
  ]
  for p in p_list:
    p.start()

  # Sum the results from each worker together
  t0 = time.time()
  xmv_mat, tot_cnt = None, 0
  for _ in range(processes):
    xmv_mat_shard, cnt = result_q.get()
    tot_cnt += cnt
    xmv_mat = (xmv_mat + xmv_mat_shard) if xmv_mat is not None else xmv_mat_shard
  t1 = time.time()
  logger.debug('Processed {} reads in {:.2f}s ({:.2f} r/s)'.format(tot_cnt, t1 - t0, tot_cnt/(t1 - t0)))

  # Orderly exit
  for p in p_list:
    p.join()

  return xmv_mat


def save(xmv_mat, fname):
  np.save(fname, xmv_mat)


def worker(worker_id, bam_fname, long_qname_table, max_xd=200, max_MQ=70, max_vlen=200, strict_scoring=False,
           in_q=None, out_q=None):
  """
  :param worker_id:
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
  logger.debug('Starting worker {} ...'.format(worker_id))

  t0, tot_cnt = time.time(), 0

  xmv_mat = np.zeros(shape=(2 * max_xd + 3, max_MQ + 1, 2 * max_vlen + 1 + 2 + 1), dtype=int)
  v_off_idx, max_v_idx = max_vlen + 2, 2 * max_vlen + 3
  bam_fp = pysam.AlignmentFile(bam_fname)
  for reference in iter(in_q.get, __process_stop_code__):
    logger.debug('Worker {}: Contig {} ...'.format(worker_id, reference))
    t1 = time.time()
    cnt = process_contig(bam_fp, long_qname_table, max_MQ, max_v_idx, max_xd, reference, strict_scoring, v_off_idx,
                         xmv_mat)
    t2 = time.time()
    logger.debug(
      'Worker {}: Contig {}: {} reads in {:2f}s ({:2f} r/s)'.format(worker_id, reference, cnt, t2 - t1, cnt / (t2 - t1)))
    tot_cnt += cnt

  out_q.put([xmv_mat, tot_cnt])

  t1 = time.time()
  logger.debug('Worker {}: Processed {} reads in {:2f}s ({:2f} r/s)'.format(worker_id, tot_cnt, t1 - t0, tot_cnt / (t1 - t0)))


def process_contig(bam_fp, long_qname_table, max_MQ, max_v_idx, max_xd, reference, strict_scoring, v_off_idx,
                   xmv_mat):
  cnt = 0
  for r in fetch(bam_fp, reference=reference):
    if r.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    cnt += 1
    ri = parse_qname(r.qname, long_qname_table=long_qname_table)[1 if r.is_read2 else 0]
    d_err = score_alignment_error(r, ri=ri, max_d=max_xd, strict=strict_scoring)
    if not ri.v_list:
      xmv_mat[max_xd + d_err, min(r.mapping_quality, max_MQ), 0] += 1
    else:
      for v_size in ri.v_list:
        xmv_mat[max_xd + d_err, min(r.mapping_quality, max_MQ), min(max_v_idx, max(1, v_off_idx + v_size))] += 1
  return cnt


def plot_figures(xmv_mat, fig_prefix=None, plot_bin_size=5):
  plot_MQ(xmv_mat, fig_prefix=fig_prefix)
  plot_AA(xmv_mat, fig_prefix=fig_prefix, plot_bin_size=plot_bin_size)
  plot_mean_MQ_heatmap(xmv_mat, fig_prefix=fig_prefix, variant_bin_size=plot_bin_size, derr_bin_size=5)


def plot_MQ(xmv_mat, fig_prefix):
  fig = plt.figure(figsize=(8, 18))
  plt.subplots_adjust(bottom=0.05, top=0.98, hspace=.3, right=0.97, left=0.18)

  ax = plt.subplot(612)
  plot_heatmap_MQ_vs_derr(ax, xmv_mat)

  ax = plt.subplot(611)
  plot_mean_MQ_vs_derr(ax, xmv_mat)

  ax = plt.subplot(613)
  plot_perr_vs_MQ(ax, xmv_mat, yscale='log')

  ax = plt.subplot(614)
  plot_perr_vs_MQ(ax, xmv_mat, yscale='linear')

  ax = plt.subplot(615)
  plot_read_count_vs_MQ(ax, xmv_mat)

  ax = plt.subplot(616)
  plot_read_fate(ax, xmv_mat)

  plt.savefig(fig_prefix + '_MQ.png')


def plot_heatmap_MQ_vs_derr(ax, xmv_mat):
  hm = xmv_mat.sum(axis=2).T
  max_derr = int((xmv_mat.shape[0] - 3) / 2)
  x1 = max_derr * 1.3
  xt = [-max_derr, -int(max_derr / 2), 0, int(max_derr / 2), max_derr]
  xlim = [-max_derr * 1.1, x1 * 1.1]

  ax.imshow(hm[:, :-2], cmap=plt.cm.gray_r, norm=LogNorm(vmin=0.01, vmax=hm.max()),
            origin='lower', extent=(xt[0], xt[-1], 0, hm.shape[0] - 1),
            aspect='auto')

  ax.set_xticks(xt)
  ax.set_xticklabels(xt)
  ax.set_xlim(xlim)

  ax.set_yticks([0, 20, 40, 60])
  ax.set_ylim([-5, 70])
  # ax.axvline(x=0, color='k', linestyle=':')

  plt.xlabel(r'$d_{err}$')
  plt.ylabel('MQ')


def plot_mean_MQ_vs_derr(ax, xmv_mat):
  """Plot mean_MQ (y-axis) against d_err (x-axis) for given range of variant sizes

  :param xmv_mat:
  :param fig_prefix:
  :param plot_bin_size:
  :return:
  """
  mq_mat = xmv_mat.sum(axis=2)
  data_cnt = mq_mat.sum(axis=1)
  mq_vector = np.arange(mq_mat.shape[1])
  mean_mq = np.dot(mq_vector, mq_mat.T) / np.clip(data_cnt, 1, data_cnt.max() + 1).astype(float)

  max_derr = int((mean_mq.shape[0] - 3) / 2)
  x1 = max_derr * 1.3
  xt = [-max_derr, -int(max_derr/2), 0, int(max_derr/2), max_derr]
  xlim = [-max_derr * 1.1, x1 * 1.1]

  ax.plot(range(-max_derr, max_derr + 1), mean_mq[:2 * max_derr + 1], 'k.')
  # ax.scatter(range(-max_derr, max_derr + 1), mean_mq[:2 * max_derr + 1],
  #            400 * data_cnt[:2 * max_derr + 1] / data_cnt[:2 * max_derr + 2].max(),
  #            color='none', edgecolors='k')
  ax.plot([x1], mean_mq[2 * max_derr + 1:2 * max_derr + 2], 'ko')
  # ax.scatter([x1], mean_mq[2 * max_derr + 1],
  #            400 * mean_mq[2 * max_derr + 1] / mean_mq[:2 * max_derr + 2].max(),
  #            color='none', edgecolors='k')

  ax.set_xticks(xt + [x1])
  ax.set_xticklabels(xt + ['WC'])
  ax.set_xlim(xlim)
  # ax.set_xlim((-50, 50))

  ax.set_yticks([0, 20, 40, 60])
  ax.set_ylim([-5, 70])
  ax.axvline(x=0, color='k', linestyle=':')

  plt.xlabel(r'$d_{err}$')
  plt.ylabel('Mean MQ')

  # plt.imshow(xmv_mat.sum(axis=2).T, cmap=plt.cm.gray_r, norm=LogNorm(vmin=0.01, vmax=1e6))


def plot_perr_vs_MQ(ax, xmv_mat, yscale='linear'):
  mq = np.arange(xmv_mat.shape[1])
  mq_mat = xmv_mat.sum(axis=2)
  for fmt, d_margin, label in zip(['g*', 'ro'], [0, 50], [r'$d_{err} = 0$', r'$|d_{err}| \leq 50$']):
    p_err = compute_p_err(mq_mat, d_margin)
    ax.plot(mq, p_err, fmt, mfc='none', label=label)

  ax.plot(mq, 10 ** (-mq / 10), 'k:')
  plt.setp(ax, xlabel='MQ', ylabel=r'$p_{err}$', xlim=(-1, 65), yscale=yscale)
  ax.legend(loc='upper right' if yscale=='linear' else 'lower left')


def compute_p_err(mq_mat, d_margin=0):
  d_max = int((mq_mat.shape[0] - 3)/2)
  c = mq_mat[d_max - d_margin:d_max + 1 + d_margin, :].sum(axis=0)
  t = mq_mat.sum(axis=0)
  return (t - c) / np.clip(t, 1, 1 + t.max())


def plot_read_count_vs_MQ(ax, xmv_mat):
  mq = np.arange(xmv_mat.shape[1])
  mq_mat = xmv_mat.sum(axis=2)
  cnt = mq_mat.sum(axis=0)
  ax.step(mq, cnt, 'k', where='mid')
  plt.setp(ax, xlabel='MQ', ylabel=r'Read count', xlim=(-1, 65), yscale='log')


def plot_read_fate(ax, xmv_mat):
  d_max = int((xmv_mat.shape[0] - 3) / 2)

  d0 = xmv_mat[d_max:d_max + 1, :, :].sum()
  within_d50 = xmv_mat[d_max - 50:d_max + 1 + 50, :, :].sum()
  within_dmax = xmv_mat[:2 * d_max + 1, :, :].sum()
  wc = xmv_mat[2 * d_max + 1, :, :].sum()
  um = xmv_mat[2 * d_max + 2, :, :].sum()

  ax.barh(0, d0, color='g', align='center')
  ax.barh(1, within_d50 - d0, color='y', align='center')
  ax.barh(2, within_dmax - within_d50, color='r', align='center')
  ax.barh(3, wc, color='pink', align='center')
  ax.barh(4, um, color='c', align='center')

  ax.text(d0, 0, d0, va='center', ha='right')
  ax.text(within_d50 - d0, 1, within_d50 - d0, va='center', ha='left')
  ax.text(within_dmax - within_d50, 2, within_dmax - within_d50, va='center', ha='left')
  ax.text(wc, 3, wc, va='center', ha='left')
  ax.text(um, 4, um, va='center', ha='left')

  ax.set_yticks([0, 1, 2, 3, 4])
  ax.set_yticklabels([r'$d_{err} = 0$', r'$0 < |d_{err}| \leq 50$', r'50 < $|d_{err}|$', 'WC', 'UM'])
  ax.set_ylim([-0.5, 4.5])
  ax.set_ylabel('Read fate')

  ax.set_xscale('log')
  ax.set_xlabel('Read count')


def plot_AA(xmv_mat, fig_prefix, plot_bin_size=5):
  fig = plt.figure(figsize=(6, 6))
  plt.subplots_adjust(bottom=0.05, top=0.98, hspace=.01)

  ax = plt.subplot(211)
  plot_alignment_accuracy_by_vsize(ax, xmv_mat, plot_bin_size=plot_bin_size)

  ax = plt.subplot(212)
  plot_vcounts(ax, xmv_mat, plot_bin_size=plot_bin_size)

  plt.savefig(fig_prefix + '_V.png')


def alignment_accuracy_by_vsize(xmv_mat, d_err_threshold):
  max_derr = int((xmv_mat.shape[0] - 3) / 2)
  derr_mat = xmv_mat.sum(axis=1)
  d_idx = (max_derr - d_err_threshold, max_derr + d_err_threshold + 1)
  num = derr_mat[d_idx[0]:d_idx[1], :].sum(axis=0)[1:]
  den = derr_mat.sum(axis=0)[1:]
  return num, den


def plot_alignment_accuracy_by_vsize(ax, xmv_mat, plot_bin_size=5):
  """

  :param ax:
  :param xmv_mat:
  :param plot_bin_size:
  :param show_ax_label:
  :return:
  """
  n0, d0 = alignment_accuracy_by_vsize(xmv_mat, 0)
  n50, d50 = alignment_accuracy_by_vsize(xmv_mat, 50)

  d0_p = plot_panels(ax,
                    num=n0,
                    den=d0,
                    bin_size=plot_bin_size,
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='b', label=r'$d_{err} = 0$')
  d50_p = plot_panels(ax,
                    num=n50,
                    den=d50,
                    bin_size=plot_bin_size,
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='r', label=r'$|d_{err}| \leq 50$')
  plt.legend(handles=[d0_p, d50_p], loc='lower center', fontsize=9)
  ax.set_ylabel('Fraction correctly aligned')


def plot_vcounts(ax, xmv_mat, plot_bin_size=5):
  v_cnt = xmv_mat.sum(axis=0).sum(axis=0)[1:]
  n_max = 10 ** np.ceil(np.log10(v_cnt.max() + 2))
  n_min = 10 ** int(np.log10(v_cnt.min() + 1))
  tot_p = plot_panels(ax,
                      num=v_cnt,
                      bin_size=plot_bin_size,
                      yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                      color='k', label='Count', show_ax_label=True, yscale='log')
  plt.legend(handles=[tot_p], loc='lower center', fontsize=9)
  ax.set_ylabel('Read count')


def plot_mean_MQ_heatmap(xmv_mat, fig_prefix, variant_bin_size=5, derr_bin_size=5):
  binned_MQ, binned_rc, v_bins, d_bins = bin_MQ_xmv_2d(xmv_mat, variant_bin_size=variant_bin_size, derr_bin_size=derr_bin_size)

  fig = plt.figure(figsize=(12, 6))
  plt.subplots_adjust(bottom=0.1, left=0.1, top=0.98, right=0.95, hspace=0.1)

  ax1 = plt.subplot(121)
  plt.imshow(binned_MQ, cmap=plt.cm.gray_r,
            origin='lower', extent=(v_bins[0], v_bins[-1], d_bins[0], d_bins[-1]),
            aspect='auto', interpolation='none')
  plt.colorbar(label='MQ')
  plt.xlabel('Variant size')
  plt.ylabel(r'$d_{err}$')

  ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)  # plt.subplot(122)
  plt.imshow(binned_rc, cmap=plt.cm.gray_r,
            origin='lower', extent=(v_bins[0], v_bins[-1], d_bins[0], d_bins[-1]),
            aspect='auto', interpolation='none',
            norm=LogNorm(vmin=0.5, vmax=binned_rc.max() + 1))
  plt.colorbar(label='Read count')
  plt.xlabel('Variant size')

  plt.savefig(fig_prefix + '_MQ_hm.png')


def bin_MQ_xmv_2d(xmv_mat, variant_bin_size=5, derr_bin_size=20):
  """Bin the xmv_mat in cells variant_bin_size x derr_bin_size and fill it with mean MQ for each cell

  :param xmv_mat:
  :param variant_bin_size:
  :param derr_bin_size:
  :return:
  """
  b_v, b_d = variant_bin_size, derr_bin_size
  d_max = int((xmv_mat.shape[0] - 3) / 2)
  n_d = max(2, int(d_max / b_d))
  v_max = int((xmv_mat.shape[2] - 3) / 2)
  n_v = max(2, int(v_max / b_v))

  mq_vector = np.arange(xmv_mat.shape[1])

  binned_MQ, binned_read_count, v_bins, d_bins = \
    np.zeros(shape=(2 * n_d + 1, 2 * n_v + 1), dtype=float), np.zeros(shape=(2 * n_d + 1, 2 * n_v + 1), dtype=float), [], []
  for i in range(-n_d, n_d + 1):
    for j in range(-n_v, n_v + 1):
      d_slice, v_slice = slice(i, b_d), slice(j, b_v)
      mq_mat = xmv_mat[d_max + d_slice[0]:d_max + d_slice[1], :, v_max + 2 + v_slice[0]:v_max + 2 + v_slice[1]].sum(axis=0).sum(axis=1)
      binned_read_count[n_d + i, n_v + j] = mq_mat.sum()
      binned_MQ[n_d + i, n_v + j] = (mq_vector * mq_mat).sum() / max(1, mq_mat.sum())

      if i == 0:
        v_bins += [(v_slice[0] + v_slice[1] - 1) / 2]
      if j == 0:
        d_bins += [(d_slice[0] + d_slice[1] - 1) / 2]

  return binned_MQ, binned_read_count, v_bins, d_bins


def slice(k, b_k):
  k_slice = (0, 1)
  if k < 0:
    k_slice = (b_k * k, b_k * (k + 1))
  elif k > 0:
    k_slice = (b_k * (k - 1) + 1, b_k * k + 1)
  return k_slice
