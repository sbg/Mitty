import os

import mitty.lib.vcfio as vio

import mitty.tests


def basic_tests():
  """I/O: Basic vcf loading."""
  for fn in ['g0_s0.vcf', 'g0_s0.vcf.gz']:
    df = vio.read_sample_from_vcf(os.path.join(mitty.tests.example_data_dir, fn), 'g0_s0')
    assert df.columns[1] == 'POS', df.columns
    assert df.ix[0].GT == 2
