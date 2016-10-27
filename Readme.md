Quickstart
==========

Run tests
---------

```
nosetests mitty.test -v
```

Data note
---------
The BED files and other small data files used in this example are found under examples/reads in the source
distribution. The example assume that the commands are being run from inside this directory. 


### VCF

The original sample VCF for the examples below is the Genome In a Bottle truth data set for 
NA12878_HG001/NISTv3.3.1/GRCh37 obtained from [the official NCBI ftp site][giab]. I assume that this
file has been saved under the name hg001.vcf.gz in the working directory. The bed file that is used
(hg001.bed) selects out 10MB portions of chromosome 1 and chromosome 10.

[giab]: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/



Prepare VCF file for taking reads from
--------------------------------------
The read simulator can not properly generate ground truth reads from complex variants (variant calls where the REF and
ALT are both larger than 1 bp). Such calls must be filtered from a VCF file before it can be used as a sample to
generate reads from.

```
mitty -v4 filter-variants hg001.vcf.gz INTEGRATION hg001.bed - 2> filter.log | bgzip -c > hg001.filt.vcf.gz
tabix -p vcf hg001.filt.vcf.gz
```

`filter.log` looks like::

```
$ cat filter.log 
DEBUG:mitty.lib.vcfio:Starting filtering ...
DEBUG:mitty.lib.vcfio:Filtering ('1', 20000, 1020000)
DEBUG:mitty.lib.vcfio:Filtered out 1:943126 CTTTTTTTTTTTTTTTTTTTTTTTT -> ('C', 'CTTTTTTTTTT')
DEBUG:mitty.lib.vcfio:Filtering ('10', 20000, 1020000)
DEBUG:mitty.lib.vcfio:Filtered out 10:169616 GAAGGA -> ('GAGGA', 'G')
DEBUG:mitty.lib.vcfio:Filtered out 10:190580 CTT -> ('C', 'CT')
DEBUG:mitty.lib.vcfio:Filtered out 10:195225 GTTATTATTA -> ('GTTATTATTATTA', 'G')
DEBUG:mitty.lib.vcfio:Filtered out 10:206889 CTTTTTTTTTTT -> ('C', 'CTTTTT')
DEBUG:mitty.lib.vcfio:Filtered out 10:276412 CAAAAA -> ('C', 'CA')
DEBUG:mitty.lib.vcfio:Filtered out 10:307381 TAAA -> ('TA', 'T')
DEBUG:mitty.lib.vcfio:Filtered out 10:396572 ATTT -> ('A', 'AT')
DEBUG:mitty.lib.vcfio:Filtered out 10:450150 ATT -> ('A', 'ATTT')
DEBUG:mitty.lib.vcfio:Filtered out 10:513505 TA -> ('T', 'TAAA')
DEBUG:mitty.lib.vcfio:Filtered out 10:572470 CAAA -> ('CAA', 'C')
DEBUG:mitty.lib.vcfio:Processed 1475 variants
DEBUG:mitty.lib.vcfio:Filtered out 11 complex variants
DEBUG:mitty.lib.vcfio:Took 0.20363688468933105 s
```

**NOTE: If the BED file is not sorted, the output file needs to be sorted again.**

The file `hg001.filt.vcf.gz` can now be used to generate reads and will serve as a truth VCF for VCF benchmarking.



Listing read models
-------------------

Listing Mitty's built in models
```
mitty list-read-models
```

Listing additional custom models stored in a folder somewhere
```
mitty list-read-models -d ./
```


Prepare a Illumina type read model from a BAM
---------------------------------------------
**This is only needed if one of the existing read models does not match your requirements. This allows you to 
sample reads from a BAM and build an empirical read model for Illumina**
_This assumes that a sample BAM file has been downloaded to the working directory_

```
mitty -v4 bam2illumina sample.bam ./rd-model.pkl "This model is taken from a BAM of unknown provenance" --every 5
```

  
Generating perfect reads
------------------------
```
mitty -v4 generate-reads ~/Data/human_g1k_v37_decoy.fasta hg001.filt.vcf.gz INTEGRATION hg001.bed ./rd-model.pkl 30 7 >(gzip > r1.fq.gz) --fastq2 >(gzip > r2.fq.gz) --threads 2
```

_(When you supply a model file name to the read generator it will first look among the builtin
read models to see if the file name is a match (typically these are in the `mitty/data/readmodels`
folder). It will then treat the model file name as a path and try and load that from your file system 
- which is the case in this particular example.)_


Corrupting reads
----------------
The reads generated using the previous command have no base call errors. Base call errors can be introduced into the
reads using the following command.

```
mitty -v4 corrupt-reads ./rd-model.pkl r1.fq.gz >(gzip > r1c.fq.gz) 7 --fastq2-in r2.fq.gz --fastq2-out >(gzip > r2c.fq.gz)
```

Note that using a read model different to that used to generate the reads originally can lead to undefined behavior, 
including crashes.


Combining read generation and corruption using named pipes
----------------------------------------------------------
```
mkfifo tf1
mkfifo tf2
mitty -v4 generate-reads ~/Data/human_g1k_v37_decoy.fasta hg001.filt.vcf.gz INTEGRATION hg001.bed ./rd-model.pkl 30 7 >(tee > tf1 >(gzip > r1.fq.gz)) --fastq2 >(tee > tf2 >(gzip > r2.fq.gz)) --threads 2 &
mitty -v4 corrupt-reads ./rd-model.pkl tf1 >(gzip > r1c.fq.gz) 7 --fastq2-in tf2 --fastq2-out >(gzip > r2c.fq.gz)
```

**Note: On Darwin reading from named pipes is a bit broken and the following construction is needed instead**
```
mitty -v4 corrupt-reads ./rd-model.pkl <(cat tf1) >(gzip > r1c.fq.gz) 7 --fastq2-in <(cat tf2) --fastq2-out >(gzip > r2c.fq.gz)
```


Alignment with BWA
------------------

```
bwa mem ~/Data/human_g1k_v37_decoy.fasta r1c.fq.gz r2c.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwac.bam
samtools index bwac.bam
```


These alignments are easy to inspect using a genome browser

![IGV screenshot showing read qname and one het variant](docs/images/igv-alignment-qname.png?raw=true "IGV screenshot showing read qname and one het variant")

Since the qname carries the correct alignment and CIGAR string you can match that against the actual alignment and
CIGAR string for spot checks.


Alignment diagnostics
---------------------

## Mapping quality
```
mitty -v4 mq-plot bwac.bam bwac.mq.csv bwac.mq.png
```


### Strict scoring
```
mitty -v4 mq-plot bwac.bam bwac.mq.csv bwac.mq.png --strict-scoring
```



## Aignment error analysis
```
mitty -v4 derr-plot bwac.bam bwac.derr.csv bwac.derr.png
```


Perfect BAM (God aligner)
-------------------------

```
mitty -v4 god-aligner ~/Data/human_g1k_v37_decoy.fasta r1c.fq.gz perfect.bam --fastq2 r2c.fq.gz --threads 4
```

Find differences in alignments
------------------------------
(This requires bamUtils to be installed. I found release 1.0.14 to work, but repository head NOT to)

```
bam diff --in1 bwac.bam --in2 perfect.bam --out diffd.unsorted.bam
samtools sort diffd.unsorted.bam > diffd.bam
samtools index diffd.bam
```
