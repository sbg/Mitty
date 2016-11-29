#!/usr/bin/env bash
set -x

FASTA=$1
VCF=$2
SAMPLE=$3
BED=$4
MODEL=$5
COVERAGE=$6
SEED=$7
FILE_PREFIX=$8
THREADS=$9

mitty -v4 generate-reads \
  ${FASTA} \
  ${VCF} \
  ${SAMPLE} \
  ${BED} \
  ${MODEL} \
  ${COVERAGE} \
  ${SEED} \
  >(gzip > ${FILE_PREFIX}_r1.fq.gz) \
  ${FILE_PREFIX}_lq.txt \
  --fastq2 >(gzip > ${FILE_PREFIX}_r2.fq.gz) \
  --threads ${THREADS}

sleep 2s  # https://kaushikghose.wordpress.com/2016/11/23/sleep-on-it/