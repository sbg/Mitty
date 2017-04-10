#!/usr/bin/env bash

set -xe

SAMPLEVCF=hg001.small.vcf.gz
SAMPLENAME=HG001
REGION_BED=region.bed
FILTVCF=hg001.filt.test.vcf

# Select sample and region from VCF file
mitty -v4 filter-variants \
  ${SAMPLEVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  ${FILTVCF} \
  2> vcf-filter.log

cat ${FILTVCF}

cat << EOF
Chrom 1 should start with a 0/0 entry followed by several normal entries, chrom 3 should have only one 0/0 entry.

This checks
1. Even though chrom 3 has two regions, we only save the 0/0 once.

* The reason why chrom 1 has a 0/0 entry is because the first region is empty
EOF