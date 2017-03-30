#!/usr/bin/env bash
# The CINS model creates insertions that are copies of the reference genome. These insertions are
# challenging to align to

set -x
mitty -v4 simulate-variants \
  cins.vcf \
  ~/Data/human_g1k_v37_decoy.fasta \
  S1 \
  region.bed \
  7 \
  --p-het 0.6 \
  --model CINS 0.0001 100 1000