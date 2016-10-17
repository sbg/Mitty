import os

from nose.tools import assert_raises, assert_sequence_equal

import mitty.lib.vcfio as vio
import mitty.test


def test_basic():
  """I/O: Basic vcf loading."""
  v = vio.load_variant_file(
    os.path.join(mitty.test.example_data_dir, 'tiny.vcf.gz'), 'g0_s0',
    os.path.join(mitty.test.example_data_dir, 'tiny.8-14.bed'))

  # A test of pysam properly filtering variants to region
  assert v[0]['v'][1][0].tuple() == (11, 'CAA', 'C', 'D', 2), v[0]['v']
  assert v[0]['v'][0][0].tuple() == (14, 'G', 'T', 'X', 0), v[0]['v']
  assert len(v[0]['v'][0]) == 1


def test_complex_variant_error():
  """I/O: Flag complex variants in VCF"""
  assert_raises(ValueError,
                vio.load_variant_file,
                os.path.join(mitty.test.example_data_dir, 'flawed-tiny.vcf.gz'),
                'g0_s0',
                os.path.join(mitty.test.example_data_dir, 'tiny.whole.bed'))


def load_data():
  """A set of reference and variants for testing read generation."""
  seq = open(os.path.join(mitty.test.example_data_dir, 'tiny.fasta')).readlines()[1]
  # df = vio.read_sample_from_vcf(os.path.join(mitty.test.example_data_dir, 'tiny.vcf'), 'g0_s0')
  df = vio.load_variant_file(
    os.path.join(mitty.test.example_data_dir, 'tiny.vcf.gz'),
    'g0_s0',
    os.path.join(mitty.test.example_data_dir, 'tiny.whole.bed'))
  return seq, df


def test_snp_expansion1():
  """I/O: SNP expansion basic"""
  ref_seq, vcf = load_data()
  v = vcf[0]['v'][1][0]
  assert v.cigarop == 'X', v
  assert v.oplen == 0, v


def test_ins_expansion1():
  """I/O: INS expansion basic"""
  ref_seq, vcf = load_data()
  v = vcf[0]['v'][1][1]
  assert v.cigarop == 'I', v
  assert v.oplen == 3, v


def test_del_expansion1():
  """I/O: DEL expansion basic"""
  ref_seq, vcf = load_data()
  v = vcf[0]['v'][1][2]
  assert v.cigarop == 'D', v
  assert v.oplen == 2, v