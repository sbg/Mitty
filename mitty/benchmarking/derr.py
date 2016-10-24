"""Go through BAM file, scoring each read and classifying it by variants carried"""
from multiprocessing import Pool
import time
import logging

import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm

from mitty.simulation.reads import parse_qname


logger = logging.getLogger(__name__)


def process_bam(bam_fname, max_v_len=200, max_d=200, processes=2, out_csv=None):
  bam_fp = pysam.AlignmentFile(bam_fname)
  if max_v_len < 51:
    logger.warning('Setting max_v_len to 51')
    max_v_len = 51
  d_err_mat = np.zeros((2 * max_v_len + 2, 2 * max_d + 1), dtype=np.uint64)
  # The last row is for reference only reads

  t_st = t0 = time.time()
  tot_rd = 0
  p = Pool(processes=processes)
  for shard in p.imap_unordered(process_bam_part_w,
                                ((bam_fname, s['SN'], max_v_len, max_d) for s in bam_fp.header['SQ'])):
    rd_cnt = shard.sum()
    tot_rd += rd_cnt
    d_err_mat += shard
    dt = time.time() - t0
    logger.debug('{} templates processed ({:0.2f} templates/sec).'.format(rd_cnt, rd_cnt / dt))
    t0 = time.time()

  t_nd = time.time()
  logger.debug('{} templates processed in {:0.2f}s ({:0.2f} templates/sec).'.format(tot_rd, t_nd - t_st, tot_rd / (t_nd - t_st)))

  if out_csv is not None:
    logger.debug('Saving d_err matrix to {}'.format(out_csv))
    np.savetxt(out_csv, d_err_mat, delimiter=',', fmt='%d')

  return d_err_mat


def process_bam_part_w(args):
  return process_bam_part(*args)


def process_bam_part(bam_fname, reference, max_v_len=200, max_d=200):
  """Work out MQ and D for

  :param bam_fname:
  :param reference:
  :param max_v_len
  :param max_d:
  :return:
  """
  d_err_mat = np.zeros((2 * max_v_len + 2, 2 * max_d + 1), dtype=np.uint64)
  bam_fp = pysam.AlignmentFile(bam_fname)
  for r in bam_fp.fetch(reference=reference):
    ri = parse_qname(r.qname)[1 if r.is_read2 else 0]  # TODO: Will this work for SE reads?
    d_err = max(min((r.pos + 1 - ri.pos if r.reference_name == ri.chrom else max_d), max_d), -max_d)
    if len(ri.v_list):
      for v in ri.v_list:
        if 0 <= max_v_len + v < 2 * max_v_len + 1:
          d_err_mat[max_v_len + v, d_err + max_d] += 1
    else:
      d_err_mat[-1, d_err + max_d] += 1  # Reference read
  return d_err_mat


def plot_derr(derr_mat, plot):
  """
  1. Plot a heat map of d_err vs event size
  2. Plot d_err distributions for Ref reads, SNPs, short, medium and long INS/DEL separately

  :param derr_mat:
  :param plot:
  :return:
  """
  v_max = int((derr_mat.shape[0] - 2) / 2)
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99, left=0.15, right=0.97, hspace=0.05)

  variant_cats = [
    ('Ref', (2 * v_max + 1, 2 * v_max + 2)),
    ('SNP', (v_max, v_max + 1)),
    ('DEL (1 - 10bp)', (v_max - 10, v_max)),
    ('DEL (11 - 50 bp', (v_max - 50, v_max - 10)),
    ('DEL (> 50 bp', (0, v_max - 50)),
    ('INS (1 - 10bp)', (v_max + 1, v_max + 11)),
    ('INS (11 - 50 bp', (v_max + 11, v_max + 51)),
    ('INS (> 50 bp)', (v_max + 51, 2 * v_max + 1)),
  ]

  for n, vcat in enumerate(variant_cats):
    ax = plt.subplot(len(variant_cats), 1, n + 1)
    plot_derr_for_v(derr_mat, vcat, ax,
                    xlabel=None if n < len(variant_cats) - 1 else r'$d_{err}$',
                    ylabel=None if n < len(variant_cats) - 1 else 'Read\n distribution')

  plt.savefig(plot)


def plot_derr_for_v(derr_mat, vcat, ax, xlabel=None, ylabel=None):
  derr_max = (derr_mat.shape[1] - 1)/2
  x = np.arange(-derr_max, derr_max + 1, dtype=int)

  rd_cnt = derr_mat[vcat[1][0]:vcat[1][1], :].sum()
  ax.plot(x, derr_mat[vcat[1][0]:vcat[1][1], :].sum(axis=0) / max(1, rd_cnt))
  plt.setp(ax, yscale='log', ylim=[1e-5, 1.0])
  if xlabel is not None:
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
  else:
    plt.setp(ax, xticklabels=[], yticklabels=[])
  ax.text(x[-1], 0.9, vcat[0], va='top', ha='right')