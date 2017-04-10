#!/usr/bin/env bash

set -xe

SAMPLEVCF=human-m-ref.vcf.gz
SAMPLENAME=ref
REGION_BED=region2.bed
FILTVCF=ref.filt.test.vcf

# Select sample and region from VCF file
mitty -v4 filter-variants \
  ${SAMPLEVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  ${FILTVCF} \
  2> vcf-filter.log

cat ${FILTVCF}

cat << EOF
There should be only one entry per contig and these should be 0/0

This checks
1. Even though chrom 3 has two regions, we only save the 0/0 once.
2. Even though the single empty variants in each chrom are outside the regions, they are included
EOF