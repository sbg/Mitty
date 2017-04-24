#!/usr/bin/env bash

set -xe

FASTA=~/Data/human_g1k_v37_decoy.fasta
FASTQ_PREFIX=hg001-reads
BAM=hg001-bwa.bam
GODBAM=hg001-god.bam

# Alignment with BWA
bwa mem \
  ${FASTA} \
  ${FASTQ_PREFIX}-corrupt1.fq.gz \
  ${FASTQ_PREFIX}-corrupt2.fq.gz | samtools view -bSho temp.bam
samtools sort temp.bam > ${BAM}
samtools index ${BAM}

# Alignment analysis

mitty -v4 debug alignment-analysis process\
  ${BAM} \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${BAM}.alignment.npy \
  --fig-prefix ${BAM}.alignment \
  --max-d 200 \
  --max-size 50 \
  --plot-bin-size 5


# God aligner
mitty -v4 god-aligner \
  ${FASTA} \
  ${FASTQ_PREFIX}-corrupt1.fq.gz \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${GODBAM} \
  --fastq2 ${FASTQ_PREFIX}-corrupt2.fq.gz \
  --threads 2


