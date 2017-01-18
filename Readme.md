Mitty is a data simulator meant to help debug aligners and variant callers. 
It is released under the [Apache](LICENSE.txt) license.

![Mitty read simulation](docs/images/mitty-reads.png?raw=true)

![Mitty read simulation](docs/images/analysis-flow.png?raw=true)

Features
========

- Generate reads from the whole genome, a single tiny region or from a set of regions as desired
- Handles X, Y chromosomes and polyploidy IF the VCF GT field is properly set
- Read qname stores correct POS, CIGAR and the sizes of variants covered by the read
- Name of sample included in read allowing us to mix FASTQs from different simulations/sources
    - Can mix viral contamination into reads
    - Can do tumor/normal mixes
- Corruption module adds sequencing errors to reads
    - Read models can be sampled from existing BAM files
- "God aligner" writes out a BAM with perfect alignments which can be used for BAM comparisons

_Mitty also supplies a VCF simulator. This is currently under development_

It is also informative to browse the [release notes](release_notes.txt)

Quickstart
==========

Install
-------
```
conda create -n mymitty python=3.5
```
* Using a virtual env is recommended, but it's a personal choice
* Mitty requires Python3 to run

For external users, please use the public repository

```
pip install git+https://github.com/sbg/Mitty.git
```


For internal users testing beta, unreleased versions, please use the internal bitbucket repository

```
pip install git+https://kghose@bitbucket.sbgenomics.com/scm/ben/mitty3.git
```

Don't forget to replace `kghose` with your username


Run tests
---------

```
nosetests mitty.test -v
```


Program help
------------
Help is available from the command line. Simply invoking the base program with no arguments
will list all the sub programs with some short help:

```
mitty
```

Detailed help on particular commands is also available:

```
mitty generate-reads --help
```

Running commands with the verbose option allows you to tune what messages (ranging from Errors to Debug) you get.

```
mitty -v{1,2,3,4} <command>
```

Example simulation script 
-------------------------
For the impatient, a short but complete script for performing read simulations is given [here](examples/reads/run.sh)


Detailed tutorial with commentary
=================================

Generating reads
----------------

### Data note

The BED files and other small data files used in this example are found under examples/reads in the source
distribution. The examples assume that the commands are being run from inside this directory and the hg19 fasta is
available in a directory called `~/Data/human_g1k_v37_decoy.fasta` which also carries the index for the fasta. 
The commands are found in the `examples/reads/run.sh`

### VCF

The original sample VCF for the examples below is the Genome In a Bottle truth data set for 
NA12878_HG001/NISTv3.3.1/GRCh37 obtained from [the official NCBI ftp site][giab]. I assume that this
file has been saved under the name hg001.vcf.gz in the working directory. The bed file that is used
(hg001.bed) selects out 1MB portions of chromosome 1 and chromosome 7.

[giab]: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/


### Prepare VCF file for taking reads from
The read simulator can not properly generate ground truth reads from complex variants (variant calls where the REF and
ALT are both larger than 1 bp) and from any entry that does not precisely lay out the ALT sequence in the ALT column,
such as angle bracket ID references and variants breakend notation. 
Such calls must be filtered from a VCF file before it can be used as a sample to generate reads from.

```
mitty -v4 filter-variants hg001.vcf.gz INTEGRATION hg001.bed - 2> sim-vcf-filter.log | bgzip -c > sim-filt.vcf.gz
tabix -p vcf sim-filt.vcf.gz
```

`filter.log` looks like::

```
$ cat sim-vcf-filter.log 
DEBUG:mitty.lib.vcfio:Starting filtering ...
DEBUG:mitty.lib.vcfio:Filtering ('1', 20000, 1020000)
DEBUG:mitty.lib.vcfio:Complex variant 1:943126 CTTTTTTTTTTTTTTTTTTTTTTTT -> ('C', 'CTTTTTTTTTT')
DEBUG:mitty.lib.vcfio:Filtering ('7', 32000000, 33000000)
DEBUG:mitty.lib.vcfio:Complex variant 7:32095932 TTCTCTCTCTCTCTCTCTC -> ('T', 'TTCTCTCTCTCTCTCTCTCTC')
DEBUG:mitty.lib.vcfio:Complex variant 7:32152848 AATAT -> ('A', 'AATATATATATATATATATATATAT')
DEBUG:mitty.lib.vcfio:Complex variant 7:32168518 AATATATATATAT -> ('AATATATATATATATATATATAT', 'A')
DEBUG:mitty.lib.vcfio:Complex variant 7:32335505 CAAA -> ('C', 'CA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32367381 CAA -> ('C', 'CA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32378622 GCACA -> ('G', 'GCA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32428769 CGT -> ('C', 'CGTGTGTGTGTGT')
DEBUG:mitty.lib.vcfio:Complex variant 7:32562762 CAAA -> ('C', 'CA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32577070 GTTTTT -> ('G', 'GTTTT')
DEBUG:mitty.lib.vcfio:Complex variant 7:32632164 CTTTT -> ('CTTTTT', 'C')
DEBUG:mitty.lib.vcfio:Complex variant 7:32648707 ATTTTTTTTTT -> ('A', 'ATTT')
DEBUG:mitty.lib.vcfio:Complex variant 7:32835051 CT -> ('C', 'CTTTTT')
DEBUG:mitty.lib.vcfio:Complex variant 7:32837427 TA -> ('T', 'TAA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32972783 CAAAAAAAAAAAAAAAAA -> ('C', 'CAAAAAAAA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32980512 CAAA -> ('C', 'CAAAAA')
DEBUG:mitty.lib.vcfio:Complex variant 7:32983713 CTTTT -> ('CTT', 'C')
DEBUG:mitty.lib.vcfio:Processed 2210 variants
DEBUG:mitty.lib.vcfio:Sample had 2210 variants
DEBUG:mitty.lib.vcfio:Discarded 17 variants
DEBUG:mitty.lib.vcfio:Took 0.2465670108795166 s
```

**NOTE: If the BED file is not sorted, the output file needs to be sorted again.**

The file `hg001.filt.vcf.gz` can now be used to generate reads and will serve as a truth VCF for VCF benchmarking.


### Listing and inspecting read models

Listing Mitty's built in models
```
mitty list-read-models
```

Listing additional custom models stored in a folder somewhere
```
mitty list-read-models -d ./
```

The model statistics can be inspected by printing a panel of plots showing read characteristics
```
mitty describe-read-model 1kg-pcr-free.pkl model.png
```
See later for a list of read models supplied with Mitty and their characteristics


### Prepare a Illumina type read model from a BAM

**This is only needed if one of the existing read models does not match your requirements. This allows you to 
sample reads from a BAM and build an empirical read model for Illumina**
_The example assumes that a sample BAM file (sample.bam) has been downloaded to the working directory_

```
mitty -v4 bam2illumina sample.bam ./rd-model.pkl "This model is taken from a BAM of unknown provenance" --every 5 --min-mq 30
```

  
### Generating perfect reads

```
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
```

_(When you supply a model file name to the read generator it will first look among the builtin
read models to see if the file name is a match (typically these are in the `mitty/data/readmodels`
folder). It will then treat the model file name as a path and try and load that from your 
file system - which is the case in this particular example.)_

The produced FASTQs have qnames encoding the correct read alignments for each read in the template. qnames may exceed
the SAM specification limit (254 characters). In such cases the qname in the FASTQ is truncated to 254 characters and
the complete qname is printed in the side-car file (`lq.txt` in the example). 

The qname format can be obtained by executing `mitty qname`


#### Reference reads
As you might expect, by passing a VCF with a sample having no variants (see `human-m-ref.vcf` or `human-f-ref.vcf` under
`examples/reads`) we can generate reads with no variants, representing the reference genome.
Note the use of `0/0` for all the autosomes and the `0` for the X and Y chromosomes in the male
to indicate the proper ploidy to the simulator via these VCF files.

```
mitty -v4 generate-reads \
  ~/Data/human_g1k_v37_decoy.fasta \
  human-m-ref.vcf.gz \
  ref \
  human-male.bed \
  1kg-pcr-free.pkl \
  0.01 \
  7 \
  >(gzip > ref-r1.gq.gz) \
  ref-lq.txt \
  --fastq2 >(gzip > ref-r2.fq.gz) \
  --threads 2
```


### Corrupting reads

The reads generated using the previous command have no base call errors. Base call errors can be introduced into the
reads using the following command.

```
mitty -v4 corrupt-reads \
  1kg-pcr-free.pkl \
  r1.fq.gz >(gzip > r1c.fq.gz) \
  lq.txt \
  lqc.txt \
  7 \
  --fastq2-in r2.fq.gz \
  --fastq2-out >(gzip > r2c.fq.gz) \
  --threads 2
```

_Using a read model different to that used to generate the reads originally can lead to undefined behavior, 
including crashes._

As mentioned, the side-car file `lq.txt` carries qnames longer than 254 characters. The output side-car file `lqc.txt`
similarly carries longer qnames from the corrupted reads. The qname for each corrupted template is identical to the 
original, uncorrupted template, except for the addition of an MD-like tag that allows recovery of the original bases
before sequencing errors were introduced. The qnames in `lq.txt` can be a subset of those in `lqc.txt`.


The BQ profile of the sample FASTQ generated by this command (generated by FASTQC) looks like

Mate 1:
![FASTQC screenshot showing BQ distribution](docs/images/1kg-pcr-free-corrupt-fastqc-r1.png?raw=true "FASTQC screenshot showing BQ distribution")

Mate 2:
![FASTQC screenshot showing BQ distribution](docs/images/1kg-pcr-free-corrupt-fastqc-r2.png?raw=true "FASTQC screenshot showing BQ distribution")

can be compared with the empirical model profile shown in the Appendix for the `1kg-pcr-free.pkl` model.


### Alignment with BWA

_(Assumes bwa and samtools are installed)_
```
bwa mem ~/Data/human_g1k_v37_decoy.fasta r1c.fq.gz r2c.fq.gz | samtools view -bSho out.bam
samtools sort out.bam > bwac.bam
samtools index bwac.bam
```


These alignments are easy to inspect using a genome browser

![IGV screenshot showing read qname and one het variant](docs/images/igv-alignment-qname.png?raw=true "IGV screenshot showing read qname and one het variant")

Since the qname carries the correct alignment and CIGAR string you can match that against the actual alignment and
CIGAR string for spot checks.


### Alignment diagnostics


#### Mapping quality
```
mitty -v4 mq-plot bwac.bam bwac.mq.csv bwac.mq.png
```

![MQ analysis](docs/images/bwac.mq.png?raw=true "MQ analysis")


##### Strict scoring
The normal alignment scoring mode is to compute the difference between the aligned read position and each breakpoint in
the correct read alignment and then take the minimum of these differences. This gives credit for a 'partial' alignment
Typically this happens when the read covers a larg(er) insertion or deletion and the aligner is able to correctly place
one anchor (side of the read), but not the other one.

If strict scoring is turned on then only the first breakpoint i.e. the perfectly correct alignment is considered.


```
mitty -v4 mq-plot bwac.bam bwac.mq.csv bwac.mq.png --strict-scoring
```

![MQ analysis](docs/images/bwac.mq-strict.png?raw=true "MQ analysis")

Note the very asymmetric nature of the mean MQ plot about d_err = 0. This is because any read that starts with a 
softclip (not due to an insertion) will have it's POS value to the right of the strictly correct alignment POS.


#### Alignment error analysis
```
mitty -v4 derr-plot bwac.bam bwac.derr.csv bwac.derr.png
```

![Alignment analysis](docs/images/bwac.derr.png?raw=true "Alignment analysis")



### Perfect BAM (God aligner)

Passing the simulated FASTQ through the god aligner produces a "perfect BAM" which can be used as a truth BAM
for comparing alignments from different aligners. This truth BAM can also be used to test variant callers by
removing one moving part (the aligner) from the analysis chain.

```
mitty -v4 god-aligner \
  ~/Data/human_g1k_v37_decoy.fasta \
  r1c.fq.gz \
  lqc.txt \
  perfectc.bam \
  --fastq2 r2c.fq.gz \
  --threads 2
```

### Find differences in alignments

(This requires bamUtils to be installed. I found release 1.0.14 to work, but repository head NOT to)

```
bam diff --in1 bwac.bam --in2 perfect.bam --out diffd.unsorted.bam
samtools sort diffd.unsorted.bam > diffd.bam
samtools index diffd.bam
```


Generating samples (genomes)
----------------------------
Mitty also has features to generate simulated genomes in the form of VCF files. 

### Sampled Genomes

The `sampled-genome` command can, given a VCF representing all variants in a population and their allele 
frequencies, return a diploid sample VCF based on random sampling of the main list.

**Populations:** 

  The program is passed an `info-af` parameter that looks for an INFO field of that name. The
  program interprets the value of that field as the alternative allele frequency and uses that to sample the
  respective variant. If you have a file, say like that from the 1000G project, that represents allele frequencies
  from multiple populations for each variant, you can simulate individuals from these different populations
  by selecting the appropriate tags. e.g. for the 1000G VCFs, `EUR_AF` for Europeans, `AMR_AF` for Americans etc.

  If this parameter is omitted the `af` parameter needs to be supplied. This sets a flat alternative allele 
  frequency for all variants that is used for sampling.


**Ploidy and BED files:** 

  The program accepts multiple BED files. A chromosome may appear zero or more times in each bed file. 
  The ploidy of each simulated chromosome is the number of BED files in which it appears at least once.

  Say, for example, we have bed files that look like

```
1.bed:
1 10  1000
2 10  5000

2.bed:
1 2000  3000

3.bed:
2 1000 2000
```

If all three bed files are passed this would result in a simulated genome with two copies of both chrom 1 and chrom 2.
Due to the regions in the BED files, chrom1 would contain only HET variants while chrom2 may carry HOM variants in the
overlap region 1000-2000

If only 1.bed and 2.bed were passed the resulting simulated genome would have two copies of chrom 1 but one copy of
chrom 2


Appendix
========

Built-in read models
--------------------
![](docs/images/1kg-pcr-free.png?raw=true)

![](docs/images/hiseq-X-v2.5-Garvan.png?raw=true)

![](docs/images/old-Garvan.png?raw=true)

![](docs/images/hiseq-2500-v1-pcr-free.png?raw=true)

![](docs/images/hiseq-X-v1-HLI.png?raw=true)