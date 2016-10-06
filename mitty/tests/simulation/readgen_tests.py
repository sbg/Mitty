"""Tests for the generation of reads by applying variants to a reference"""
import os

from nose.tools import assert_raises

import mitty.tests
import mitty.simulation.readgen as rgen
import mitty.lib.vcfio as vio


def v_type_test():
  """Read gen: Flag complex variants in VCF"""
  df = vio.read_sample_from_vcf(os.path.join(mitty.tests.example_data_dir, 'flawed-tiny.vcf'), 'g0_s0')
  for _, v in df.iterrows():
    assert_raises(RuntimeError, rgen.get_v_type, v)


def load_data():
  """A set of reference and variants for testing read generation."""
  seq = open(os.path.join(mitty.tests.example_data_dir, 'tiny.fasta')).readlines()[1]
  df = vio.read_sample_from_vcf(os.path.join(mitty.tests.example_data_dir, 'tiny.vcf'), 'g0_s0')
  return seq, df


def snp_expansion_test1():
  """Read gen: SNP expansion basic"""
  ref_seq, df = load_data()
  v = df.ix[0]
  vt, vl = rgen.get_v_type(v)
  assert vt == 'X', vt
  assert vl == 0, vl


def snp_expansion_test2():
  """Read gen: SNP expansion: No M section"""
  ref_seq, df = load_data()
  v = df.ix[0]
  samp_pos = 1
  ref_pos = 5
  vt, vl = rgen.get_v_type(v)
  nodes, n_samp_pos, n_ref_pos = rgen.snp(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 1, nodes
  assert nodes[0] == (1, 5, 'X', 1, 'T', 0), nodes[0]
  assert n_samp_pos == samp_pos + 1
  assert n_ref_pos == ref_pos + 1


def snp_expansion_test3():
  """Read gen: SNP expansion: M section"""
  ref_seq, df = load_data()
  v = df.ix[0]
  samp_pos = 1
  ref_pos = 1
  vt, vl = rgen.get_v_type(v)
  nodes, samp_pos, ref_pos = rgen.snp(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 2, nodes
  assert nodes[0] == (1, 1, '=', 4, 'ATGA', None), nodes
  assert nodes[1] == (5, 5, 'X', 1, 'T', 0), nodes


def ins_expansion_test1():
  """Read gen: INS expansion basic"""
  ref_seq, df = load_data()
  v = df.ix[1]
  vt, vl = rgen.get_v_type(v)
  assert vt == 'I', vt
  assert vl == 3, vl


def ins_expansion_test2():
  """Read gen: INS expansion: No M section"""
  ref_seq, df = load_data()
  v = df.ix[1]
  samp_pos = 9
  ref_pos = 9
  vt, vl = rgen.get_v_type(v)
  nodes, n_samp_pos, n_ref_pos = rgen.insertion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 1, nodes
  assert nodes[0] == (9, 9, 'I', 3, 'TTT', 3), nodes[0]
  assert n_samp_pos == samp_pos + 3
  assert n_ref_pos == ref_pos


def ins_expansion_test3():
  """Read gen: INS expansion: M section"""
  ref_seq, df = load_data()
  v = df.ix[1]
  samp_pos = 6
  ref_pos = 6
  vt, vl = rgen.get_v_type(v)
  nodes, samp_pos, ref_pos = rgen.insertion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 2, nodes
  assert nodes[0] == (6, 6, '=', 3, 'GTA', None), nodes[0]
  assert nodes[1] == (9, 9, 'I', 3, 'TTT', 3), nodes[1]


def del_expansion_test1():
  """Read gen: DEL expansion basic"""
  ref_seq, df = load_data()
  v = df.ix[2]
  vt, vl = rgen.get_v_type(v)
  assert vt == 'D', vt
  assert vl == 2, vl


def del_expansion_test2():
  """Read gen: DEL expansion: No M section"""
  ref_seq, df = load_data()
  v = df.ix[2]
  samp_pos = 12
  ref_pos = 12
  vt, vl = rgen.get_v_type(v)
  nodes, n_samp_pos, n_ref_pos = rgen.deletion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 1, nodes
  assert nodes[0] == (12, 14, 'D', 2, '', -2), nodes[0]
  assert n_samp_pos == samp_pos, n_samp_pos
  assert n_ref_pos == ref_pos + 2, n_ref_pos


def del_expansion_test3():
  """Read gen: DEL expansion: M section"""
  ref_seq, df = load_data()
  v = df.ix[2]
  samp_pos = 12
  ref_pos = 9
  vt, vl = rgen.get_v_type(v)
  nodes, samp_pos, ref_pos = rgen.deletion(ref_seq, samp_pos, ref_pos, v, vl)

  assert len(nodes) == 2, nodes
  assert nodes[0] == (12, 9, '=', 3, 'TCC', None), nodes
  assert nodes[1] == (15, 14, 'D', 2, '', -2), nodes


def expand_sequence_test():
  """Sequence expansion"""
  ref_seq, df = load_data()
  nodes = rgen.create_node_list(ref_seq, ref_start_pos=1, chrom_copy=0b01, vcf=df)

  assert len(nodes) == 9
  #                   ps pr op l seq
  assert nodes[0] == (1, 1, '=', 4, 'ATGA', None)
  assert nodes[1] == (5, 5, 'X', 1, 'T', 0)
  assert nodes[2] == (6, 6, '=', 3, 'GTA', None)
  assert nodes[3] == (9, 9, 'I', 3, 'TTT', 3)
  assert nodes[4] == (12, 9, '=', 3, 'TCC', None)
  assert nodes[5] == (15, 14, 'D', 2, '', -2)
  assert nodes[6] == (15, 14, '=', 7, 'GGAGGCG', None)
  assert nodes[7] == (22, 23, 'D', 2, '', -2)
  assert nodes[8] == (22, 23, '=', 3, 'ACC', None)
