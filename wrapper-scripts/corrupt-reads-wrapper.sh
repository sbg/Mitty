#!/usr/bin/env bash
set -x

MODEL=$1
R1=$2
LQ=$3
R2=$4
SEED=$5
FILE_PREFIX=$6
THREADS=$7

mitty -v4 corrupt-reads \
  ${MODEL} \
  ${R1} \
  >(gzip > ${FILE_PREFIX}_r1c.fq.gz) \
  ${LQ} \
  ${FILE_PREFIX}_lqc.txt \
  ${SEED} \
  --fastq2-in ${R2} \
  --fastq2-out >(gzip > ${FILE_PREFIX}_r2c.fq.gz) \
  --threads ${THREADS}

sleep 2s  # https://kaushikghose.wordpress.com/2016/11/23/sleep-on-it/