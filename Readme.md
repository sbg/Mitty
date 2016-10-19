Quickstart
==========

Run tests
---------

```
nosetests mitty.test -v
```

Prepare VCF file for taking reads from
--------------------------------------
(Run from inside examples/hg19)

Filter out complex variants from file, restrict to required sample and given BED file

```
mitty -v4 filter-variants NA12878.vcf.gz NA12878 hg19_1_Y_breakpoints.bed - 2> filter.log | bgzip -c > NA12878.filtered.vcf.gz
tabix -p vcf NA12878.filtered.vcf.gz
```

`filter.log` looks like::

```
  $ tail filter.log 
  DEBUG:mitty.lib.vcfio:Filtered out 22:50423000 GT -> ('GTT', 'GTT')
  DEBUG:mitty.lib.vcfio:Filtered out 22:50899755 AGGG -> ('AGG', 'AGG')
  DEBUG:mitty.lib.vcfio:Filtered out 22:51028063 CAAACA -> ('CAAACAAAACA', 'CAAACAAAACA')
  DEBUG:mitty.lib.vcfio:Filtered out 22:51108762 CA -> ('CAA', 'CAA')
  DEBUG:mitty.lib.vcfio:Filtering ('X', 60001, 115682290)
  DEBUG:mitty.lib.vcfio:Filtering ('X', 115732291, 155270560)
  DEBUG:mitty.lib.vcfio:Filtering ('Y', 2649521, 59373566)
  DEBUG:mitty.lib.vcfio:Processed 4101558 variants
  DEBUG:mitty.lib.vcfio:Filtered out 15016 complex variants
  DEBUG:mitty.lib.vcfio:Took 58.96592593193054 s
```

**NOTE: If the BED file is not sorted, the output file needs to be sorted again.**


Taking reads
------------

The Illumina parameters are stored in a `.json file` and look like:

```
  {
    "gc_bias": null,
    "rlen": 150,
    "tlen": 500,
    "tlen_std": 30,
    "diploid_coverage": 1
  }
```
  
Taking reads looks like::

```
mitty -v4 generate-reads illumina150x500.json ~/Data/human_g1k_v37_decoy.fasta NA12878.filtered.vcf.gz NA12878 hg19_10_breakpoints.bed 7 --fastq1 >(gzip > r1.fq.gz) --fastq2 >(gzip > r2.fq.gz) --threads 2
```

Alignment with BWA
------------------

```
bwa mem ~/Data/human_g1k_v37_decoy.fasta r1.fq.gz r2.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > sorted.bam
samtools index sorted.bam
```

Perfect BAM (God aligner)
-------------------------

```
mitty -v4 god-aligner ~/Data/human_g1k_v37_decoy.fasta r1.fq.gz perfect.bam --fastq2 r2.fq.gz --threads 4
```



Empirical Base Quality Score
----------------------------

```
gunzip -c ~/Downloads/sim-reads.fq.gz | mitty -v4 sample-fastq-base-quality --max-reads 100000 -t 2 - - > test.npz
```