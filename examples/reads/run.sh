#!/usr/bin/env bash
set -ex


#mitty -v4 bam2illumina sample.bam ./rd-model.pkl "This model is taken from a BAM of unknown provenance" --every 5 --min-mq 30
#mitty describe-read-model ./rd-model.pkl model.png


mitty -v4 filter-variants hg001.vcf.gz INTEGRATION region.bed - 2> filter.log | bgzip -c > hg001.filt.vcf.gz
tabix -p vcf hg001.filt.vcf.gz

mitty -v4 generate-reads \
  ~/Data/human_g1k_v37_decoy.fasta \
   hg001.filt.vcf.gz \
   INTEGRATION \
   region.bed \
   1kg-pcr-free.pkl \
   10 \
   7 \
   >(gzip > r1.fq.gz) \
   lq.txt \
   --fastq2 >(gzip > r2.fq.gz) \
   --threads 2

mitty -v4 corrupt-reads \
  1kg-pcr-free.pkl \
  r1.fq.gz >(gzip > r1c.fq.gz) \
  lq.txt lqc.txt \
  7 \
  --fastq2-in r2.fq.gz \
  --fastq2-out >(gzip > r2c.fq.gz) \
  --threads 2

bwa mem ~/Data/human_g1k_v37_decoy.fasta r1c.fq.gz r2c.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwac.bam
samtools index bwac.bam

mitty -v4 debug alignment-analysis-plot bwac.bam lqc.txt xmv-data

exit 0

mitty -v4 god-aligner \
  ~/Data/human_g1k_v37_decoy.fasta \
  r1.fq.gz \
  lq.txt \
  perfect.bam \
  --fastq2 r2.fq.gz \
  --platform-name Illumina \
  --sample-name INTEGRATION \
  --threads 2

mitty -v4 god-aligner \
  ~/Data/human_g1k_v37_decoy.fasta \
  r1c.fq.gz \
  lqc.txt \
  perfectc.bam \
  --fastq2 r2c.fq.gz \
  --sample-name INTEGRATION \
  --threads 2


bam diff --in1 bwac.bam --in2 perfect.bam --out diffd.unsorted.bam
samtools sort diffd.unsorted.bam > diffd.bam
samtools index diffd.bam

# Cleanup
# rm *.fq.gz *.bam *.bai *.png *.log
