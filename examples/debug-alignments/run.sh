#!/usr/bin/env bash
set -ex
# Generate a set of reads (FASTQ) and then run BWA MEM three times on the same data with different parameters
# Then compare the BAMs read by read using the bam-partition tool

mitty -v4 generate-reads \
  ~/Data/human_g1k_v37_decoy.fasta \
   hg001.filt.vcf.gz \
   INTEGRATION \
   hg001.bed \
   1kg-pcr-free.pkl \
   30 \
   7 \
   >(gzip > r1.fq.gz) \
   lq.txt \
   --fastq2 >(gzip > r2.fq.gz) \
   --threads 2

# Generate three BAMs from the same reads by running BWA MEM with different parameters
l=("1.5" "10" "20")
for i in ${l[*]}
do
  bwa mem -r ${i} ~/Data/human_g1k_v37_decoy.fasta r1.fq.gz r2.fq.gz | samtools view -bSho out.bam
  samtools sort out.bam > bwa_${i}.bam
  samtools index bwa_${i}.bam
done

bwa mem -k 50 ~/Data/human_g1k_v37_decoy.fasta r1.fq.gz r2.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwa_seed50.bam
samtools index bwa_seed50.bam


bwa mem -w 1 ~/Data/human_g1k_v37_decoy.fasta r1.fq.gz r2.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwa_w1.bam
samtools index bwa_w1.bam



# Run the BAM partitioner
mitty -v4 debug partition-bams \
  derr \
  d_err --threshold 10 \
  --sidecar_in lq.txt --bam bwa_1.5.bam --bam bwa_10.bam --bam bwa_20.bam

mitty -v4 debug partition-bams \
  map mapped \
  --bam bwa_1.5.bam  --bam bwa_10.bam --bam bwa_20.bam

mitty -v4 debug partition-bams \
  mq MQ \
  --threshold 1  --bam bwa_1.5.bam --bam bwa_10.bam --bam bwa_20.bam