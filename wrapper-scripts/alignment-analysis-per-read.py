"""
This script takes an indexed BAM file created from simulated reads and analyzes
them for alignment accuracy.

Commandline is of the form

  python3 alignment-analysis-per-read.py  <BAM> <LONG-QNAME> <OUT-DIR> <NCPUS>


The output is a tab delimited file with the following fields:

chrom pos d_err MQ  variantlist

  variantlist is a semicolon separated, spaceless list of variant sizes that looks like:
  -1;+2;+1...

The output should be sorted, bgzipped and then tabix indexed for maximum
usefulness, such as extracting particular bam regions etc.

An example workflow is

python3 alignment-analysis-per-read.py  alignments.bam  long-qname.txt ./ 4 \
| sort -k 1,1 -k 2,2 -n -t $'\t' \
| bgzip > scored-reads.gz \
&& tabix -s 1 -b 2 -e 2 scored-reads.gz

"""
import sys
import os
from subprocess import call

import mitty.analysis as maly


bam_fname = sys.argv[1]
scar_fname = sys.argv[2]
output_dir = sys.argv[3]
ncpus = int(sys.argv[4])


max_d = 200

n1 = maly.parse_read_qnames(scar_fname)
n2 = maly.compute_derr(max_d=max_d)
n3 = maly.save_as_tab(maly.get_seq_dict(bam_fname), output_dir)

pipeline = [n1, n2, n3]


output_file_list = list(maly.scatter(pipeline, bam_fname=bam_fname, ncpus=ncpus))

script = "cat " + " ".join(output_file_list)
call(script, shell=True)

for fn in output_file_list:
  os.remove(fn)