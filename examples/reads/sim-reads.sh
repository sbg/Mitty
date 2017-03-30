#!/usr/bin/env bash
set -xe

FASTA=~/Data/human_g1k_v37_decoy.fasta
SAMPLEVCF=hg001.vcf.gz
SAMPLENAME=HG001
REGION_BED=region.bed
FILTVCF=hg001-filt.vcf.gz
READMODEL=1kg-pcr-free.pkl
COVERAGE=30
READ_GEN_SEED=7
FASTQ_PREFIX=hg001-reads
READ_CORRUPT_SEED=8

# Select sample and region from VCF file
mitty -v4 filter-variants \
  ${SAMPLEVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  - \
  2> vcf-filter.log | bgzip -c > ${FILTVCF}

tabix -p vcf ${FILTVCF}

# Take reads
mitty -v4 generate-reads \
  ${FASTA} \
  ${FILTVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  ${READMODEL} \
  ${COVERAGE} \
  ${READ_GEN_SEED} \
  >(gzip > ${FASTQ_PREFIX}1.fq.gz) \
  ${FASTQ_PREFIX}-lq.txt \
  --fastq2 >(gzip > ${FASTQ_PREFIX}2.fq.gz) \
  --threads 2

# Corrupt reads
mitty -v4 corrupt-reads \
  ${READMODEL} \
  ${FASTQ_PREFIX}1.fq.gz >(gzip > ${FASTQ_PREFIX}-corrupt1.fq.gz) \
  ${FASTQ_PREFIX}-lq.txt \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${READ_CORRUPT_SEED} \
  --fastq2-in ${FASTQ_PREFIX}2.fq.gz \
  --fastq2-out >(gzip > ${FASTQ_PREFIX}-corrupt2.fq.gz) \
  --threads 2