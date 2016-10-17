Cookbook
========

Run tests
---------

nosetests mitty.test -v

Prepare VCF file for taking reads from
--------------------------------------

::

  mitty -v4 filter-variants ../examples/hg19/NA12878.vcf.gz NA12878 ../examples/hg19/hg19_1_Y_breakpoints.bed - | bgzip -c > ../examples/hg19/NA12878.filtered.vcf.gz
  tabix -p vcf NA12878.filtered.vcf.gz

Output looks like::

  DEBUG:mitty.lib.vcfio:Filtered out 22:49391407 TCACA -> ('TCACACACA', 'TCACA')
  DEBUG:mitty.lib.vcfio:Filtered out 22:49638783 AATGGATGGATGG -> ('AATGG', 'AATGG')
  DEBUG:mitty.lib.vcfio:Filtered out 22:50290460 AT -> ('AT', 'ATT')
  DEBUG:mitty.lib.vcfio:Filtered out 22:50423000 GT -> ('GTT', 'GTT')
  DEBUG:mitty.lib.vcfio:Filtered out 22:50899755 AGGG -> ('AGG', 'AGG')
  DEBUG:mitty.lib.vcfio:Filtered out 22:51028063 CAAACA -> ('CAAACAAAACA', 'CAAACAAAACA')
  DEBUG:mitty.lib.vcfio:Filtered out 22:51108762 CA -> ('CAA', 'CAA')
  DEBUG:mitty.lib.vcfio:Filtered out 15016 complex variants

**NOTE: If the BED file is not sorted, the output file needs to be sorted again.**


mitty -v4 filter-variants ../examples/hg19/HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_CHROM1-22_v3.2.2_highconf.vcf.gz INTEGRATION ../examples/hg19/hg19_1_Y_breakpoints.bed - | bgzip -c > ../examples/hg19/HG002.filtered.vcf.gz

Taking reads
------------

mitty -v4 generate-reads illumina150x500.json ~/Data/human_g1k_v37_decoy.fasta ../hg19/HG002.filtered.vcf.gz INTEGRATION ../hg19/hg19_10_breakpoints.bed 7 --fastq1 r1.fq --fastq2 r2.fq --threads 2


mitty -v4 generate-reads illumina150x500.json ~/Data/human_g1k_v37_decoy.fasta ../hg19/HG002.filtered.vcf.gz INTEGRATION ../hg19/hg19_10_breakpoints.bed 7 --fastq1 >(gzip > r1.fq.gz) --fastq2 >(gzip > r2.fq.gz) --threads 2


Empirical Base Quality Score
----------------------------

gunzip -c ~/Downloads/sim-reads.fq.gz | mitty -v4 sample-fastq-base-quality --max-reads 100000 -t 2 - - > test.npz