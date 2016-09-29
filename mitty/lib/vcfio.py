import gzip

import numpy as np
import pandas as pd


_gt_dict = {
  '0|0': 0,
  '0|1': 1,
  '1|0': 2,
  '1|1': 3,
  '0/1': 4,  # We will alternate copy for thsi
  '1/1': 2
}

def _convert_gt(gt_string):
  """Function that converts the gt string (0|1, 1|0, 1|1) to a
  two bit genotype code. If the data is het un-phased it is
  alternately assigned a chromosome.

  :param row:
  :return:
  """
  gt = _gt_dict.get(gt_string[:3])  # We assume the GT value is the first entry in FORMAT
  if gt == 4:
    gt = _convert_gt.copy + 1
    _convert_gt.copy = not _convert_gt.copy
  return gt


def read_sample_from_vcf(fname, sample):
  gz = fname.endswith('gz')
  my_open = gzip.open if gz else open
  with my_open(fname, 'r') as f:
    for n, line in enumerate(f):
      if gz:
        line = line.decode()  # gzip returns bytes
      if line.startswith('#CHROM'):
        if sample not in line:
          raise RuntimeError('Sample {} does not exist in {}'.format(sample, fname))
        break
    else:
      raise RuntimeError('Malformed VCF {}'.format(fname))

  # TODO: handle large CSV files in chunks?
  return pd.read_csv(
    fname, delimiter='\t', header=n,
    usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample],
    converters={sample: _convert_gt},
    dtype={'POS': np.uint32}).rename(columns={'#CHROM': 'CHROM', sample: 'GT'})
