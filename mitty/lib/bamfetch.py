"""Current pysam does now allow us to fetch the unmapped reads tucked away at the end of the BAM.
We fake this by creating a custom fetch generator which, if passed the last contig in the BAM first reads
those alignments, and then uses until_eof=True to read rest of the reads in the file, which are
(the unmapped reads all collected just past that)"""
from itertools import chain


def fetch(bam_fp, reference):
  return chain(*[bam_fp.fetch(reference=reference)] +
                ([bam_fp.fetch(until_eof=True)] if reference == bam_fp.references[-1] else []))
