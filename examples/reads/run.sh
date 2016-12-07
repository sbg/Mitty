#!/usr/bin/env bash
set -x

PREFIX=sim
FASTA=~/Data/human_g1k_v37_decoy.fasta
RAWVCF=hg001.vcf.gz
FILTVCF=${PREFIX}-filt.vcf.gz
SAMPLE=INTEGRATION
BED=hg001.bed
RMODEL=hiseq-2500-v1-pcr-free.pkl
COV=30
RSEED=7
CSEED=8
R1=${PREFIX}-r1.fq.gz
LQ=${PREFIX}-lq.txt
R2=${PREFIX}-r2.fq.gz
RC1=${PREFIX}-r1c.fq.gz
LQC=${PREFIX}-lqc.txt
RC2=${PREFIX}-r2c.fq.gz
BWABAM=${PREFIX}-bwa.bam
PERBAM=${PREFIX}-perfect.bam
DIFFBAM=${PREFIX}-diff.bam

#mitty -v4 bam2illumina sample.bam ./rd-model.pkl "This model is taken from a BAM of unknown provenance" --every 5 --min-mq 30
#mitty describe-read-model ./rd-model.pkl model.png

mitty -v4 filter-variants ${RAWVCF} ${SAMPLE} ${BED} - 2> ${PREFIX}-vcf-filter.log | bgzip -c > ${FILTVCF}
tabix -p vcf ${FILTVCF}

mitty -v4 generate-reads \
  ${FASTA} \
  ${FILTVCF} \
  ${SAMPLE} \
  ${BED} \
  ${RMODEL} \
  ${COV} \
  ${RSEED} \
  >(gzip > ${R1}) \
  ${LQ} \
  --fastq2 >(gzip > ${R2}) \
  --threads 2

mitty -v4 corrupt-reads \
  ${RMODEL} \
  ${R1} \
  >(gzip > ${RC1}) \
  ${LQ} \
  ${LQC} \
  ${CSEED} \
  --fastq2-in ${R2} \
  --fastq2-out >(gzip > ${RC2}) \
  --threads 2

bwa mem \
  ${FASTA} \
  ${RC1} \
  ${RC2} | samtools view -bSho out.bam
samtools sort out.bam > ${BWABAM}
samtools index ${BWABAM}


mitty -v4 god-aligner \
  ${FASTA} ${RC1} ${LQC} ${PERBAM} --fastq2 ${RC2} --threads 2

bam diff --noCigar --in1 ${BWABAM} --in2 ${PERBAM} --out out.bam
samtools sort out.bam > ${DIFFBAM}
samtools index ${DIFFBAM}

rm out.bam

# Cleanup
# rm ${PREFIX}*