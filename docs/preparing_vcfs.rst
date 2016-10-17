Preparing VCFs for the read simulator
+++++++++++++++++++++++++++++++++++++

The read simulator consumes standard VCFs but some requirements must be met if things are to work as expected:


1. If the VCF being fed to the read generator will also be used as a truth VCF then it has to be cut to the same
   regions as the BED file passed to the generator. This can be easily done by using the ``mitty filter-variants``
   command.
2. If any regions of the genome within the BED file have no variants you are required to place a dummy variant there
   with the correct GT field so that Mitty can infer the ploidy of the region correctly. If this is omitted,
   Mitty will default to generating diploid reads for that region.


Examples of dummy variants indicating ploidy are::

  ##fileformat=VCFv4.1
  ##contig=<ID=1,length=23>
  ##contig=<ID=2,length=23>
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	g0_s0
  1	5	.	C	T	100	PASS	.	GT	0|0
  2	8	.	A	ATTT	100	PASS	.	GT	0


Examples:

1. Reads from a human female are being simulated:
   The BED file should not include chrom Y. It is not necessary to have dummy variants on any chromosome since the
   default is to assume diploidy. An empty VCF file still needs to be passed.
2. Reference reads from a human male are being simulated:
   The BED file should include chrom Y. If no variants are present on Y, or if purely reference reads are being taken
   a dummy variant with GT 0 should be placed for both X and Y, otherwise coverage for X and Y will be double that
   expected


Ploid BED
---------

This is similar to a regular BED data, except that we enhance it to be able to represent individual copies of the
chromosome. Regions in a ploid BED have an additional "ploid field" indicating which chromosome copy the region
refers to. This is used internally by the simulator and allows simulating polyploid organisms and the X and Y
chromosomes for females and males properly.


