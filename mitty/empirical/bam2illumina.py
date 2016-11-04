"""Given a BAM file extract

1. Template len and std model
2. BQ model
3. GC bias model"""
import time
from multiprocessing import Pool
import pickle

import numpy as np
import pysam

import logging
logger = logging.getLogger(__name__)


def process_bam_section(bam_fname, reference, every=1, min_mq=0, max_bq=94, max_bp=300, max_tlen=1000):
  """

  :param bam_fname:
  :param reference:
  :param every:
  :param min_mq:
  :param max_bq:
  :param max_bp:
  :param max_tlen:
  :return:
  """
  bq_mat = np.zeros((2, max_bp, max_bq), dtype=np.uint64)
  tlen_mat = np.zeros(max_tlen, dtype=np.uint64)
  t0 = time.time()
  n, min_rlen, max_rlen, mean_rlen, r_cnt = 0, 1e13, 0, 0, 0
  for n, r in enumerate(pysam.AlignmentFile(bam_fname, 'rb').fetch(reference=reference)):
    # We process every read (except secondary alignments) for BQ
    # We filter by minimum MQ for template length computations
    if n % every != 0: continue
    if r.flag > 255: continue  # Ignore secondary/split alignments

    r_cnt += 1
    bq_r = r.query_qualities[::-1] if r.flag & 0x10 else r.query_qualities  # Take care of reverse complement
    bq_mat[0 if r.is_read1 else 1, range(len(bq_r)), bq_r] += 1

    if min_rlen > r.rlen: min_rlen = r.rlen
    if max_rlen < r.rlen: max_rlen = r.rlen

    if r.mapq >= min_mq and r.mapq != 255 and r.is_read2 and abs(r.tlen) < max_tlen:
      tlen_mat[abs(r.tlen)] += 1

    if n % 100000 == 0: logger.debug('Processed {} reads'.format(n))

  t1 = time.time()
  logger.debug('Took {}s to parse {} reads ({} r/s)'.format(t1 - t0, n, n/(t1 - t0)))
  return {
    'r_cnt': r_cnt,
    'bq_mat': bq_mat,
    'tlen_mat': tlen_mat,
    'min_rlen': min_rlen,
    'max_rlen': max_rlen,
  }


def process_bam_section_w(kwargs):
  return process_bam_section(**kwargs)


def process_bam_parallel(bam_fname, pkl, model_description='Test model', every=1, min_mq=0,
                         threads=4, max_bq=94, max_bp=300, max_tlen=1000):
  """

  :param bam_fname:
  :param pkl:
  :param model_description:
  :param every:
  :param min_mq:
  :param threads:
  :param max_bq:
  :param max_bp:
  :param max_tlen:
  :return:
  """
  p = Pool(threads)
  t0 = time.time()
  bam_fp = pysam.AlignmentFile(bam_fname, mode='rb')
  refs = [sq['SN'] for sq in bam_fp.header['SQ']][:24]

  bq_mat = np.zeros((2, max_bp, max_bq), dtype=np.uint64)
  tlen_mat = np.zeros(max_tlen, dtype=np.uint64)
  min_rlen, max_rlen, rlen_sum, r_cnt = 1e13, 0, 0, 0
  for data in p.imap_unordered(
    process_bam_section_w,
    ({"bam_fname": bam_fname,
      "reference": ref,
      "every": every,
      "min_mq": min_mq,
      "max_bq": max_bq,
      "max_bp": max_bp,
      "max_tlen": max_tlen}
     for ref in refs)):
    bq_mat += data['bq_mat']
    tlen_mat += data['tlen_mat']

    r_cnt += int(data['r_cnt'])

    min_rlen = min(min_rlen, data['min_rlen'])
    max_rlen = max(max_rlen, data['max_rlen'])

    t1 = time.time()
    logger.debug('Processed {} templates ({} t/s)'.format(r_cnt, r_cnt / (t1 - t0)))

  cum_bq_mat = np.zeros(bq_mat.shape, dtype=float)
  for p in [0, 1]:
    cum_bq_mat[p, :, :] = bq_mat[p, :, :].cumsum(axis=1) / bq_mat[p, :, :].sum(axis=1).clip(1)[:, None]  # clip ensures no division by zero nonsense
  # This array can be used to generate a BQ value from a uniformly distributed random value r
  # np.searchsorted(n_bq[base_no, :], r), where base_no is the base of the read being considered

  pickle.dump({
    'model_class': 'illumina',
    'model_description': model_description,
    'min_mq': min_mq,
    'bq_mat': bq_mat,
    'cum_bq_mat': cum_bq_mat,
    'tlen': tlen_mat,
    'r_cnt': r_cnt,
    'cum_tlen': tlen_mat.cumsum() / tlen_mat.sum(),  # This array can similarly be used as np.searchsorted(ct, r)
    'mean_rlen': max_rlen,  # This value is needed for coverage calculations
                            # and for this model, is also the read length used
    'min_rlen': min_rlen,
    'max_rlen': max_rlen
  }, open(pkl, 'wb'))