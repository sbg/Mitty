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

Don't forget to replace 'kghose' with your username


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
(hg001.bed) selects out 1MB portions of chromosome 1 and chromosome 10.

[giab]: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/

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



Listing and inspecting read models
----------------------------------

Listing Mitty's built in models
```
mitty list-read-models
```

Listing additional custom models stored in a folder somewhere
```
mitty list-read-models -d ./
```

The model statistics can be inspected
```
mitty describe-read-model ./rd-model.pkl model.png
```

![Example of model statistics](docs/images/model.png?raw=true "Example of model statistics")


Prepare a Illumina type read model from a BAM
---------------------------------------------
**This is only needed if one of the existing read models does not match your requirements. This allows you to 
sample reads from a BAM and build an empirical read model for Illumina**
_This assumes that a sample BAM file has been downloaded to the working directory_

```
mitty -v4 bam2illumina sample.bam ./rd-model.pkl "This model is taken from a BAM of unknown provenance" --every 5 --min-mq 30
```

  
Generating perfect reads
------------------------
```
mitty -v4 generate-reads ~/Data/human_g1k_v37_decoy.fasta hg001.filt.vcf.gz INTEGRATION hg001.bed ./rd-model.pkl 30 7 >(gzip > r1.fq.gz) --fastq2 >(gzip > r2.fq.gz) --threads 2
```

_(When you supply a model file name to the read generator it will first look among the builtin
read models to see if the file name is a match (typically these are in the `mitty/data/readmodels`
folder). It will then treat the model file name as a path and try and load that from your 
file system - which is the case in this particular example.)_


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


Alignment diagnostics
---------------------

## Mapping quality
```
mitty -v4 mq-plot bwac.bam bwac.mq.csv bwac.mq.png
```

![MQ analysis](docs/images/bwac.mq.png?raw=true "MQ analysis")


### Strict scoring
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


## Alignment error analysis
```
mitty -v4 derr-plot bwac.bam bwac.derr.csv bwac.derr.png
```

![Alignment analysis](docs/images/bwac.derr.png?raw=true "Alignment analysis")



Perfect BAM (God aligner)
-------------------------
Passing the simulated FASTQ through the god aligner produces a "perfect BAM" which can be used as a truth BAM
for comparing alignments from different aligners. This truth BAM can also be used to test variant callers by
removing one moving part (the aligner).

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
