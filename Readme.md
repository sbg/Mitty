Mitty is a data simulator meant to help debug aligners and variant callers. 
It is released under the [Apache](LICENSE.txt) license.

![Genome simulation](docs/images/genome-simulation.png?raw=true)

![Read simulation](docs/images/read-simulation.png?raw=true)

![Report generation](docs/images/reports.png?raw=true)

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
- Simple genome simulator to generate VCFs with SNPs and different sizes of Insertions and Deletions for aligner/caller testing

It is also informative to browse the [release notes](release_notes.txt) and poke around the detailed documentation under the [/docs](docs/) folder

Quickstart
==========

Install
-------
```
conda create -n mymitty python=3.5
```
* Using a virtual env is recommended, but it's a personal choice
* Mitty requires Python3 to run

Install Mitty from the public github repository

```
pip install git+https://github.com/sbg/Mitty.git
```


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


Detailed tutorial with commentary
=================================

Note on the files used here
---------------------------
A script containing all the commands in this tutorial is found under [`examples/tutorial/script.sh`](examples/tutorial/script.sh). Some parts of the script assume external tools such as an aligner, a variant caller and the GA4GH benchmarking tools are installed. 
I assume that the GRC37 (or HG19) FASTA and it's index is available and aliased to `$FASTA`. On my development machine, for example, `FASTA=~/Data/human_g1k_v37_decoy.fasta`.
The VCF refered to as `hg001.vcf.gz` has been downloaded from the [NIH NCBI website](ftp://ftp-trace.ncbi.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/)

Generating reads
----------------

### Aliases used
```
FASTA=~/Data/human_g1k_v37_decoy.fasta
SAMPLEVCF=hg001.vcf.gz
SAMPLENAME=HG001
REGION_BED=region.bed
FILTVCF=hg001-filt.vcf.gz
READMODEL=1kg-pcr-free.pkl
COVERAGE=30
READ_GEN_SEED=7
FASTQ_PREFIX=hg001-reads
READ_CORRUPT_SEED=8
```

### Prepare VCF file for taking reads from
The read simulator can not properly generate ground truth reads from overlapping variants 
(e.g. a deletion that overlaps a SNP on the same chromosome copy),  complex variants (variant calls where the REF and
ALT are both larger than 1 bp) or any entry that does not precisely lay out the ALT sequence in the ALT column,
such as angle bracket ID references and variants breakend notation. 
Such calls must be filtered from a VCF file before it can be used as a sample to generate reads from.

```
mitty -v4 filter-variants \
  ${SAMPLEVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  - \
  2> vcf-filter.log | bgzip -c > ${FILTVCF}

tabix -p vcf ${FILTVCF}
```

The `region.bed` file is used to select out just the parts of the VCF we would like to generate reads from. The input VCF can be a population VCF (say from the 1000G project). The `${SAMPLENAME}` indicates which sample to extract from the VCF.


`vcf-filter.log` looks like::

```
$ cat vcf-filter.log 
DEBUG:mitty.lib.vcfio:Starting filtering ...
DEBUG:mitty.lib.vcfio:Filtering ('1', 1000000, 2000000)
DEBUG:mitty.lib.vcfio:Complex variant 1:1827835 CTT -> ('CT', 'C')
DEBUG:mitty.lib.vcfio:Filtering ('3', 60830384, 61830384)
DEBUG:mitty.lib.vcfio:Complex variant 3:60835995 TTCTCTCTCTCTCTCTCTCTCTC -> ('TTCTCTCTCTCTCTCTCTCTCTCTCTCTC', 'T')
DEBUG:mitty.lib.vcfio:Complex variant 3:60897726 TAC -> ('TACAC', 'T')
DEBUG:mitty.lib.vcfio:Complex variant 3:60970457 GGTGTGT -> ('GGT', 'G')
DEBUG:mitty.lib.vcfio:Complex variant 3:61205001 ATG -> ('ATGTG', 'A')
DEBUG:mitty.lib.vcfio:Complex variant 3:61309628 GGTGTGTGTGTGTGTGT -> ('GGTGTGTGTGTGTGTGTGTGT', 'G')
DEBUG:mitty.lib.vcfio:Complex variant 3:61360782 AGTGTGTGTGT -> ('AGTGTGTGTGTGT', 'A')
DEBUG:mitty.lib.vcfio:Complex variant 3:61469726 CAAAAAAAAAAAAAAA -> ('CA', 'C')
DEBUG:mitty.lib.vcfio:Complex variant 3:61488707 TTGTG -> ('TTG', 'T')
DEBUG:mitty.lib.vcfio:Complex variant 3:61509647 TCACACACACA -> ('TCACACACACACACA', 'T')
DEBUG:mitty.lib.vcfio:Complex variant 3:61522251 TACACAC -> ('TACAC', 'T')
DEBUG:mitty.lib.vcfio:Complex variant 3:61541525 CA -> ('CAAA', 'C')
DEBUG:mitty.lib.vcfio:Complex variant 3:61633465 CAAAAAAA -> ('CAAAAAAAAAAAAA', 'C')
DEBUG:mitty.lib.vcfio:Complex variant 3:61718383 AAAATAAAT -> ('AAAAT', 'A')
DEBUG:mitty.lib.vcfio:Complex variant 3:61731724 CTTT -> ('CTT', 'C')
DEBUG:mitty.lib.vcfio:Processed 2807 variants
DEBUG:mitty.lib.vcfio:Sample had 2807 variants
DEBUG:mitty.lib.vcfio:Discarded 15 variants
DEBUG:mitty.lib.vcfio:Took 0.28023505210876465 s
```

**NOTE: If the BED file is not sorted, the output file needs to be sorted again.**

`${FILTVCF}` can now be used to generate reads and will serve as a truth VCF for VCF benchmarking.


### Listing and inspecting read models

`mitty list-read-models` will list built in read models 

`mitty list-read-models -d ./mydir` will list additional custom models stored in `mydir`

`mitty describe-read-model 1kg-pcr-free.pkl model.png` prints a panel of plots showing read characteristics

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
  ${FASTA} \
  ${FILTVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  ${READMODEL} \
  ${COVERAGE} \
  ${READ_GEN_SEED} \
  >(gzip > ${FASTQ_PREFIX}1.fq.gz) \
   ${FASTQ_PREFIX}-lq.txt \
   --fastq2 >(gzip > ${FASTQ_PREFIX}2.fq.gz) \
   --threads 2
```

_(When you supply a model file name to the read generator it will first look among the builtin
read models to see if the file name is a match (typically these are in the `mitty/data/readmodels`
folder). It will then treat the model file name as a path and try and load that from your 
file system - which is the case in this particular example.)_

The produced FASTQs have qnames encoding the correct read alignments for each read in the template. qnames may exceed the SAM specification limit (254 characters). In such cases the qname in the FASTQ is truncated to 254 characters and the complete qname is printed in the side-car file `${FASTQ_PREFIX}-lq.txt`. 

The qname format can be obtained by executing `mitty qname`


#### Reference reads
As you might expect, by passing a VCF with a sample having no variants (see `human-m-ref.vcf` or `human-f-ref.vcf` under
`examples/empty-vcfs`) we can generate reads with no variants, representing the reference genome.
Note the use of `0/0` for all the autosomes and the `0` for the X and Y chromosomes in the male
to indicate the proper ploidy to the simulator via these VCF files.


#### Truncating reads
For some experiments you might want to generate custom sized reads. `generate-reads` allows you to do 
this with the `--truncate-to` argument

```
FASTQ_PREFIX=hg001-truncated-reads
mitty -v4 generate-reads \
  ${FASTA} \
  ${FILTVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  ${READMODEL} \
  ${COVERAGE} \
  ${READ_GEN_SEED} \
  >(gzip > ${FASTQ_PREFIX}1.fq.gz) \
  ${FASTQ_PREFIX}-lq.txt \
  --fastq2 >(gzip > ${FASTQ_PREFIX}2.fq.gz) \
  --truncate-to 60 \
  --threads 2
```

This generates the same kind of reads as before, but all the reads are 60bp long, instead of their usual length. 
Naturally, you can not make reads longer than what the model originally specifies.


### Corrupting reads

The reads generated using the previous command have no base call errors. Base call errors can be introduced into the reads using the following command.

```
mitty -v4 corrupt-reads \
  ${READMODEL} \
  ${FASTQ_PREFIX}1.fq.gz >(gzip > ${FASTQ_PREFIX}-corrupt1.fq.gz) \
  ${FASTQ_PREFIX}-lq.txt \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${READ_CORRUPT_SEED} \
  --fastq2-in ${FASTQ_PREFIX}2.fq.gz \
  --fastq2-out >(gzip > ${FASTQ_PREFIX}-corrupt2.fq.gz) \
  --threads 2
```

_Using a read model different to that used to generate the reads originally can lead to undefined behavior, 
including crashes._

As mentioned, the side-car file `${FASTQ_PREFIX}-lq.txt` carries qnames longer than 254 characters. The output side-car file `${FASTQ_PREFIX}-corrupt-lq.txt` similarly carries longer qnames from the corrupted reads. The qname for each corrupted template is identical to the  original, uncorrupted template, except for the addition of an MD-like tag that allows recovery of the original bases before sequencing errors were introduced.

The BQ profile of the sample FASTQ generated by this command (generated by FASTQC) looks like

Mate 1:
![FASTQC screenshot showing BQ distribution](docs/images/1kg-pcr-free-corrupt-fastqc-r1.png?raw=true "FASTQC screenshot showing BQ distribution")

Mate 2:
![FASTQC screenshot showing BQ distribution](docs/images/1kg-pcr-free-corrupt-fastqc-r2.png?raw=true "FASTQC screenshot showing BQ distribution")

can be compared with the empirical model profile shown in the Appendix for the `1kg-pcr-free.pkl` model.


### Alignment with BWA

_(Assumes bwa and samtools are installed)_
```
FASTA=~/Data/human_g1k_v37_decoy.fasta
FASTQ_PREFIX=hg001-reads
BAM=hg001-bwa.bam

bwa mem \
  ${FASTA} \
  ${FASTQ_PREFIX}-corrupt1.fq.gz \
  ${FASTQ_PREFIX}-corrupt2.fq.gz | samtools view -bSho temp.bam
samtools sort temp.bam > ${BAM}
samtools index ${BAM}
```


These alignments are easy to inspect using a genome browser

![IGV screenshot showing read qname and one het variant](docs/images/igv-alignment-qname.png?raw=true "IGV screenshot showing read qname and one het variant")

Since the qname carries the correct alignment and CIGAR string you can match that against the actual alignment and
CIGAR string for spot checks.


### Perfect BAM (God aligner)

Passing the simulated FASTQ through the god aligner produces a "perfect BAM" which can be used as a truth BAM
for comparing alignments from different aligners. This truth BAM can also be used to test variant callers by
removing one moving part (the aligner) from the analysis chain.

```
GODBAM=hg001-god.bam

mitty -v4 god-aligner \
  ${FASTA} \
  ${FASTQ_PREFIX}-corrupt1.fq.gz \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${GODBAM} \
  --fastq2 ${FASTQ_PREFIX}-corrupt2.fq.gz \
  --threads 2
```


Analysis
--------
Mitty supplies some tools to help with benchmarking and debugging of aligner/caller pipelines.

## Alignment accuracy

```
mitty -v4 debug alignment-analysis process\
  ${BAM} \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${BAM}.alignment.npy \
  --fig-prefix ${BAM}.alignment \
  --max-d 200 \
  --max-size 50 \
  --plot-bin-size 10
```

This invocation will process `${BAM}` and summarize the alignment performance in a numpy data file (`${BAM}.alignment.npy`). It will plot them in a set of figures named with the prefix
`${BAM}.alignment`. The alignment error will be assessed upto a maximum of 200bp. The program will check variants from 50bp deletions to 50bp insertions, putting them into 10bp size bins. SNPs are always counted and placed in their own spearate bin. 

![MQ plots](docs/images/aligner-report-example-1.png?raw=true "MQ plots")

![Alignment accuracy plots](docs/images/aligner-report-example-2.png?raw=true "Alignment accuracy plots")


## Variant calling accuracy, parametrized by variant size

We can use a set of tools developed by the GA4GH consortium to compare a VCF produced by a pipeline with a truth VCF.
One of the outputs of the comparator tools is a VCF (called an evaluation VCF) where each call is annotated with by whether it is a TP, FN, FP or GT error.

`variant-call-analysis` is a program that summarizes the data in such an evaluation VCF in terms of variant size.

```
mitty -v4 debug variant-call-analysis process \
  eval.vcf.gz \
  caller.analysis.csv \
  --max-size 200 \
  --fig-file demo.png \
  --plot-bin-size 10  \
  --title 'Example'
```

This invocation will process `evcf.vcf.gz`, write the results as a comma separated file (`caller.analysis.csv`) and also plot them
in `demo.png`. 

![P/R plots](docs/images/caller-report-example.png?raw=true "P/R plots")


## Find differences in alignments

(This requires bamUtils to be installed. I found release 1.0.14 to work, but repository head NOT to)

```
bam diff --in1 bwac.bam --in2 perfect.bam --out diffd.unsorted.bam
samtools sort diffd.unsorted.bam > diffd.bam
samtools index diffd.bam
```


## Set differences of two or more BAM files derived from the same FASTQ(s)

Often, we would like to tweak alignment algorithm parameters, or even test new algorithms. Sometimes we just look at
the global effect of these tweaks, for example, by looking at variant calling performance. Often, however, we want to
examine more directly the effects of these experiments. One way to do this is to look in detail at the differences in
alignments for the same set of reads run through the different pipelines.

The Mitty `partition-bams` subtool allows us to do this. For a full example please see the script `examples/debug-alignments/run.sh`
In brief, say we have three different BAM files obtained by passing the same FASTQ through three different aligner tools.
In the example we run bwa mem with three different values of the `-f` parameter which trades off speed for accuracy.

We can then check (since this is a simulated FASTQ) how accurate the alignments were compared to each other by running

```
mitty -v4 debug partition-bams \
  myderr \
  d_err --threshold 10 \
  --sidecar_in lq.txt --bam bwa_1.5.bam --bam bwa_10.bam --bam bwa_20.bam
```

This command line asks the tool to use |d_err| < 10 as the set membership function. We are passing it three BAM files
(the file names refer to the `-r` values we passed `bwa mem` (1.5, 10 and 20)) and `lq.txt` is the sidecar file carrying
the qnames > 254 characters (as described previously). 

This tool produces a summary file `myderr_summary.txt` that looks like:

```
(A)(B)(C) 22331
(A)(B)C   234
(A)B(C)   0
(A)BC     3
A(B)(C)   0
A(B)C     0
AB(C)     208
ABC       199126
```

In this nomenclature A is the set and (A) is the complement of this set. The set labels A, B, C ... (upto a maximum of 10)
refer to the BAM files in sequence, in this case 1.5, 10 and 20. 

Thus, ABC means all the reads which have a |d_err| < 10 in all the three files. AB(C) means all the reads which have 
a |d_err| < 10 in A and B but not C, and so on. A reader familiar with Venn diagrams is referred to the chart below
for a translation of the three dimensional case to a three way Venn diagram. Higher dimensions are harder to visualize
as Venn diagrams.

![Sets to Venn diagram](docs/images/sets.png?raw=true "Sets to Venn diagram")


The tool also produces a set of files following the naming convention:

```
myderr_(A)(B)(C)_A.bam
myderr_(A)(B)(C)_B.bam
myderr_(A)(B)(C)_C.bam
myderr_(A)(B)C_A.bam
myderr_(A)(B)C_B.bam
myderr_(A)(B)C_C.bam
...
```
The first part of the name follows the convention outlined above. The trailing A, B, C refer to the orginal source BAM of
the reads. So `myderr_(A)(B)(C)_B.bam` carries reads from bam B that have |d_err| >= 10 in all the three BAMs.

An example of throwing these files up on a genome browser and inspecting them is given below

![IGV Bam Partitions](docs/images/igv-sets.png?raw=true "BAM partitions on IGV")

The criteria the `partition-bam` tool can be run on can be obtained by passing it the `--criteria` option.


Generating samples (genomes)
----------------------------
Mitty also has features to generate simulated genomes in the form of VCF files. 

### Simulated variants

The `simulate-variants` command generates a VCF with simulated variants. 
The program carries three basic models for variant simulation - SNPs, insertions and deletions and is invoked as follows:

```
mitty -v4 simulate-variants \
  - \                                      # Write the VCF to std out
  ~/Data/human_g1k_v37_decoy.fasta \       # reference
  mysample \                               # The name of the sample to add to
  region.bed \                             # region over which to generate variants
  7 \                                      # random number generator seed
  --p-het 0.6 \                            # probability for heterozygous variants
  --model SNP 0.001 1 1 \                  #  <model type> <p> <min-size> <max-size>
  --model INS 0.0001 10 100 \
  --model DEL 0.0001 10 100 | bgzip -c > sim.vcf.gz
  
tabix -p vcf sim.vcf.gz
```  
p is the probability of a variant being placed on any given base 
min-size and max-size indicate the size ranges of the variants produced. These are ignored for SNP

This VCF should be run through the `filter-variants` program as usual before taking reads. This is especially
important because the simulation can produce illegaly overlapping variants which will be taken out by this step.
  
```  
mitty -v4 filter-variants sim.vcf.gz mysample region.bed - 2> sim-vcf-filter.log | bgzip -c > sim-filt.vcf.gz
tabix -p vcf sim-filt.vcf.gz  
```

Please see `examples/variants/run.sh` for an example script.

The `CINS` model works just like the `INS` model except the insertion sequences, instead of being 
novel DNA sequences created with a Markov chain generator, are exact copies of random parts of 
the input reference genome. This creates insertions that are more challenging to align to and 
assemble, especially when their lengths start to exceed the template size of the sequencing 
technology used.


```
mitty -v4 simulate-variants \
  cins.vcf \
  ~/Data/human_g1k_v37_decoy.fasta \
  S1 \
  region.bed \
  7 \
  --p-het 0.6 \
  --model CINS 0.0001 100 1000
```



Miscellaneous utilities
-----------------------

Plot variant size distribution in VCF file:
```
mitty -v4 debug variant-by-size hg001.vcf.gz hg001.variant.size.csv --max-size 100 --title "HG001" --fig-file hg001.variant.png --plot-bin-size 5
```


Appendix
========

Qname format
------------
Read alignment and simulation metadata are stored in the qname in the following format. 

```
      @index|sn|chrom:copy|strand:pos:cigar:v1,v2,...:MD|strand:pos:cigar:v1,v2,...:MD|

Each section is separated by a  bar (`|`), sub-sections are separated by a colon (`:`) and 
their meanings  are as follows: 

index:  unique code for each template. If you must know, the index is simply a monotonic counter
        in base-36 notation
sn:     sample name. Useful if simulation is an ad-mixture of sample + contaminants or multiple samples
chrom:  chromosome id of chromosome the read is taken from
copy:   copy of chromosome the read is taken from (0, 1, 2, ...) depends on ploidy
strand: 0: forward strand, 1: reverse strand
pos:    Position of first (left-most) reference matching base in read (Same definition as for POS in BAM)
        One based
        For reads coming from completely inside an insertion this is, however, the POS for the insertion
cigar:  CIGAR string for correct alignment.
        For reads coming from completely inside an insertion, however, the CIGAR string is:
        '>p+nI' where:
           '>' is the unique key that indicates a read inside a long insertion
           'p' is how many bases into the insertion branch the read starts
           'n' is simply the length of the read
v1,v2,..: Comma separated list of variant sizes covered bly this read
            0 = SNP
            + = INS
            - = DEL
MD:     Read corruption MD tag. This is empty for perfect reads, but is filled with an MD formatted
        string for corrupted reads. The MD string is referenced to the original, perfect read, not
        a reference sequence.
```

Notes:

- The alignment information is repeated for every read in the template.
  The example shows what a PE template would look like.
- The pos value is are one based to make comparing qname info in genome browser easier
- qnames longer than 254 characters (the maximum allowed by the SAM spec) are stored in
  a side-car file alongside the simulated FASTQs.
  The qname in the FASTQ file is truncated to 254 characters. A truncated qname is detected by the
  absence of the final '|' character in the qname

This information is also available by executing `mitty qname`

Example from a perfect read:

`@8|INTEGRATION|1:1|1:1930067:27=1X63=1I158=:0,1:|0:1929910:184=1X63=1I1=:0,1:|`

Corresponding read after passing through read corruption model:

`@8|INTEGRATION|1:1|1:1930067:27=1X63=1I158=:0,1:126T29T16G15G28C0T2T8C9T0C6T|0:1929910:184=1X63=1I1=:0,1:9G25G3T1C36T87G10G0C0A4A1G4C0A0G0G0A1A0T0C3A1C0A0T0C0C0T0G2T0A0A1A0G1G0T0G0A0A0A0C1T1A0T0C0T0C0T0A1T1A0A1A0T0A0C0A0A|`


```
@8           - simulated reads serial number (in base 26)
INTEGRATION  - name of sample (from VCF) reads are generated from
1:1          - read is from chromosome 1, copy 1 (copy numbers start 0, 1, ...)
1:1930067:27=1X63=1I158=:0,1:126T29T16G15G28C0T2T8C9T0C6T 
             - read one metadata with fields as:

    1:1930067       - read is from reverse strand (1, as opposed to 0 - see the metdata for the mate) at pos 1930067
    27=1X63=1I158=  - this is the correct CIGAR for this read
    0,1             - this read carries two variants a SNP (0) and an insertion (1)
    126T29T16G15G28C0T2T8C9T0C6T  
                    - This string follows the conventions for an MD tag and indicates the sequencing error
                      relative to the uncorrupted read

0:1929910:184=1X63=1I1=:0,1:9G25G3T1C36T87G10G0C0A4A1G4C0A0G0G0A1A0T0C3A1C0A0T0C0C0T0G2T0A0A1A0G1G0T0G0A0A0A0C1T1A0T0C0T0C0T0A1T1A0A1A0T0A0C0A0A 
             - read two metadata in same format
```


Built-in read models
--------------------
![](docs/images/1kg-pcr-free.png?raw=true)

![](docs/images/hiseq-X-v2.5-Garvan.png?raw=true)

![](docs/images/old-Garvan.png?raw=true)

![](docs/images/hiseq-2500-v1-pcr-free.png?raw=true)

![](docs/images/hiseq-X-v1-HLI.png?raw=true)