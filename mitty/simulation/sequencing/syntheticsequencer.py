"""A fully synthetic read model that allows us to produce single end or paired end reads with arbitrary
read and template lengths. It's read model format is as follows

{
  'model_class': 'illumina',
  'model_description': '',
  'paired': True/False,
  'read_length': 100,
  'mean_template_length': 300,
  'std_template_length': 100,
  'bq_mat': [],
  'cum_bq_mat': []
}

"""
import pickle

import numpy as np


def create_model(
  pkl,
  read_length=100, mean_template_length=500, std_template_length=100, max_tlen=1000,
  bq0=30, k=200, sigma=10,
  comment=''):
  description = """This is a synthetic read model that generates reads
  with a length of {} bp, a template length of {} +/- {} bp.
  The mean base quality follows the equation:
        {} * exp(- b/l * {})
        where b is the base in the read and l is the length of the read.
  The base quality for a given base in a given read is drawn from a gaussian with standard deviation {}
  {}""".format(
    read_length,
    mean_template_length,
    std_template_length,
    bq0, k, sigma,
    comment)

  bq = bq0 * (1 - np.exp(- k * np.linspace(1, 0, read_length) ** 2))
  one_bq_mat = np.zeros((read_length, 94), dtype=float)
  for n in range(read_length):
    one_bq_mat[n, :] = np.exp(- 0.5 * ((np.arange(94) - bq[n]) / sigma) ** 2)
  one_cum_bq_mat = one_bq_mat.cumsum(axis=1) / one_bq_mat.sum(axis=1).clip(1)[:, None]

  tlen_mat = np.exp(- 0.5 * ((np.arange(max_tlen) - mean_template_length) / std_template_length) ** 2)
  tlen_mat /= tlen_mat.sum()
  cum_tlen = tlen_mat.cumsum() / tlen_mat.sum()

  pickle.dump({
    'model_class': 'illumina',
    'model_description': description,
    'min_mq': 0,
    'bq_mat': np.array((one_bq_mat, one_bq_mat)),
    'cum_bq_mat': np.array((one_cum_bq_mat, one_cum_bq_mat)),
    'tlen': tlen_mat,
    'cum_tlen': cum_tlen,
    'mean_rlen': read_length,
    'min_rlen': read_length,
    'max_rlen': read_length,
    'r_cnt': 1
  }, open(pkl, 'wb'))
