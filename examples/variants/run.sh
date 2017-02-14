#!/usr/bin/env bash
set -x
mitty -v4 simulate-variants \
  - \
  ~/Data/human_g1k_v37_decoy.fasta \
  S1 \
  region.bed \
  7 \
  --p-het 0.6 \
  --model SNP 0.001 1 1 \
  --model INS 0.0001 10 100 \
  --model DEL 0.0001 10 100 | bgzip -c > sim.vcf.gz

tabix -p vcf sim.vcf.gz

mitty -v4 filter-variants sim.vcf.gz S1 region.bed - 2> sim-vcf-filter.log | bgzip -c > sim-filt.vcf.gz
tabix -p vcf sim-filt.vcf.gz

mitty -v4 generate-reads \
  ~/Data/human_g1k_v37_decoy.fasta \
   sim-filt.vcf.gz \
   S1 \
   region.bed \
   1kg-pcr-free.pkl \
   30 \
   7 \
   >(gzip > r1.fq.gz) \
   lq.txt \
   --fastq2 >(gzip > r2.fq.gz) \
   --threads 2

bwa mem ~/Data/human_g1k_v37_decoy.fasta r1.fq.gz r2.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwa.bam
samtools index bwa.bam
