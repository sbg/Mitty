"""Tests for the generation of a random individual from given vcf"""
import os

import pysam
import mitty.test
import mitty.simulation.genome.sampledgenome as sg


def read_data(outname):
  vcf_pointer = pysam.VariantFile(filename=outname)
  sample = vcf_pointer.header.samples[0]
  gts = []
  for rec in vcf_pointer.fetch():
    gts.append("|".join([str(ind) for ind in rec.samples["HG"].allele_indices]))
  return sample, gts


def test_assign_random_gt():
  tempfilename = "tempfile.vcf"
  sg.assign_random_gt(input_vcf=os.path.join(mitty.test.example_data_dir, 'gt_free.vcf'),
                      output=open(tempfilename, "w"), sample_name="HG", default_af=0.1,
                      seed=5081973)

  sample_name, gt_array = read_data(tempfilename)
  os.remove(tempfilename)

  assert sample_name == "HG"
  assert "0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|1_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|1" \
         "_1|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_1|0" \
         "_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|1_0|0_0|0_0|0_0|0_0|0_0|0_1|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0" \
         "_0|0_1|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0_0|0" == "_".join(gt_array)
