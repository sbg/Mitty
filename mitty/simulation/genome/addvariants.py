"""Two common procedures we employ when simulating genomes are:
 1. Adding newly simulated variants onto an existing list representing a chromosome
    - This requires us to make sure illegal variant combinations don't occur, such as
      deletions overlapping other features
 2. "Zipping up" one or more copies of a chromosome
    - This requires resolving homozygous or polyzygous variants

This library contains routines that implement these procedures.
"""


def zip_up(vl1, vl2):
  """Given two lists of variants zip them up to produce a polyploid representation

  :param vl1: represents a genome - might have GT information
  :param vl2: represents a copy of a chromosome, GT information is ignored
  :return: new_vl1: A list of variants with ploidy = ploidy of vl1 + 1
  """

