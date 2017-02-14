# Follows standard BED convention: https://genome.ucsc.edu/FAQ/FAQformat#format1
# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# chromStart - The starting position of the feature in the chromosome or scaffold.
#              The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold.
#            The chromEnd base is not included in the display of the feature.
#
# For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100,
# and span the bases numbered 0-99.

def read_bed(bed_fname):
  return list(map(lambda x: (x[0], int(x[1]), int(x[2])), map(lambda x: x.split(), open(bed_fname, 'r').readlines())))