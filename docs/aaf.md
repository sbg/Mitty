AAF
===

The Alignment Analysis Format file is a tab delimited file with the following fields:

```
chrom pos d_err MQ  variantlist
```

`variantlist` is a semicolon separated, spaceless list of variant sizes that 
looks like `-1;+2;+1...`

This file can be indexed with tabix and read with pysam. It carries a per read
analysis of alignment and its advantages are

1. Tiny size: even though it carries a per-read analysis it is small. 
   A 1.2 GB BAM file when processed results in a 60MB `.aaf.bgz` file,
   A 113.7 GB BAM file results in a 2.1 GB aaf file
2. You can use all the normal `tabix` tools on it after indexing, a major one
   being to stratify the file using a bed file.

While scoring alignment analysis on the BAM (see `bamtoolz` and `bamfilters`) 
is very powerful in that we can extract the original reads based on different 
filters processing BAM files can be fairly unweildy. The AAF is a lightweight
subsitute that can be useful when looking at summary statistics in different
ways. After this exploratory analysis we can go back to detailed read analysis
on the original BAMs.