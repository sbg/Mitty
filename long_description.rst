Mitty is a data simulator meant to help debug aligners and variant
callers.

It requires Python 3.4 or later. It is released under the
Apache 2.0 license.

Features
========

-  Generate reads from the whole genome, a single tiny region or from a
   set of regions as desired
-  Handles X, Y chromosomes and polyploidy IF the VCF GT field is
   properly set
-  Read qname stores correct POS, CIGAR and the sizes of variants
   covered by the read
-  Name of sample included in read allowing us to mix FASTQs from
   different simulations/sources

   -  Can mix viral contamination into reads
   -  Can do tumor/normal mixes

-  Corruption module adds sequencing errors to reads

   -  Read models can be sampled from existing BAM files

-  "God aligner" writes out a BAM with perfect alignments which can be
   used for BAM comparisons
-  Simple genome simulator to generate VCFs with SNPs and different
   sizes of Insertions and Deletions for aligner/caller testing
-  Analyze, plot and debug alignment accuracy with composable analysis functions


Tutorial and usage manual can be found here_.

.. _here: https://github.com/sbg/Mitty