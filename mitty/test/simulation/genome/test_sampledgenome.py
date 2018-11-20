"""Tests for the generation of a random individual from given vcf"""
import os

import pysam
import mitty.test
import mitty.simulation.genome.sampledgenome as sg


def read_data(outname):
  vcf_pointer = pysam.VariantFile(filename=outname)
  sample = vcf_pointer.header.samples[0]
  gts = {}
  for rec in vcf_pointer.fetch():
    gts[rec.pos] = "|".join([(str(ind) if ind is not None else ".") for ind in rec.samples["HG"].allele_indices])
  return sample, gts


def test_assign_random_gt():
  """Simulation: sampled-genome"""
  tempfilename = "tempfile.vcf"

  #gt for overlapping deletion is 0|1
  sg.assign_random_gt(input_vcf=os.path.join(mitty.test.example_data_dir, 'test-sample-genome-with-af.vcf'),
                      output=open(tempfilename, "w"), sample_name="HG", default_af=0.1,
                      seed=1)

  sample_name, gt_map = read_data(tempfilename)

  assert sample_name == "HG"
  assert gt_map[30923] == "0|." or gt_map[30923] == "1|."

  # gt for overlapping deletion is 0|0
  sg.assign_random_gt(input_vcf=os.path.join(mitty.test.example_data_dir, 'test-sample-genome-with-af.vcf'),
                      output=open(tempfilename, "w"), sample_name="HG", default_af=0.1,
                      seed=2)

  _, gt_map = read_data(tempfilename)
  assert gt_map[30923] != ".|."

  # gt for overlapping deletion is 1|1
  sg.assign_random_gt(input_vcf=os.path.join(mitty.test.example_data_dir, 'test-sample-genome-with-af.vcf'),
                      output=open(tempfilename, "w"), sample_name="HG", default_af=0.1,
                      seed=204)

  _, gt_map = read_data(tempfilename)
  assert gt_map[30923] == ".|."

  # gt for overlapping deletion is 1|0
  sg.assign_random_gt(input_vcf=os.path.join(mitty.test.example_data_dir, 'test-sample-genome-with-af.vcf'),
                      output=open(tempfilename, "w"), sample_name="HG", default_af=0.1,
                      seed=124)

  _, gt_map = read_data(tempfilename)
  assert gt_map[30923] == ".|0" or gt_map[30923] == ".|1"

  os.remove(tempfilename)

  assert "0|0_0|1_0|0_0|0_0|0_0|0_0|0_0|0_0|0_1|0_0|0_0|0_0|0_1|0_0|0_0|0_0|0_0|1_0|0_0|1_0|0_0|0_0|0_0|0_0|0_0|0_0|0_" \
         "1|0_.|._1|1_0|0_0|0_.|._0|0_0|0_0|1_0|0_0|0_0|0_1|0_0|0_0|0_0|0_0|0_0|0_0|1_.|._0|0_0|0_0|0_0|0_1|0_0|0_0|0_" \
         "0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|1_0|0_0|0_0|0_0|0_0|0_0|1_0|0_0|0_1|0_.|0_0|0_0|1_0|0_1|0_1|0_0|0_0|0_0|0_0|0_" \
         "0|1_0|0_0|0_0|0_0|0_1|1_.|._0|0_1|0_0|0" == "_".join(gt_map.values())
