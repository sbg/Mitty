"""A simple routine to load in a FASTQ file and give us the distribution of qname lengths, because I was curious"""
import numpy as np
import pysam

from mitty.benchmarking.alignmentscore import load_qname_sidecar


def main(fastq_fname, qname_overflow_fname, max_expected_qname_length=500):
  long_qname_table = load_qname_sidecar(qname_overflow_fname)
  qname_count = [0] * (max_expected_qname_length + 1)
  with pysam.FastxFile(fastq_fname) as fh:
    for r in fh:
      qlen = len(long_qname_table.get(r.name.split('|', 1)[0])) if r.name[-1] != '*' else len(r.name)
      qname_count[min(qlen, max_expected_qname_length)] += 1

  return qname_count