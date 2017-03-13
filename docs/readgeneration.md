Read Generation Algorithms
==========================

Overview
--------
The algorithms are designed to load up only the parts of the reference and sample variants designated by a
supplied BED file. This allows us scale simulations smoothly in terms of memory and compute time requirements.
For example, for quick tests, we can take high coverage data from a small part of the genome, using modest
memory and time requirements. By simply switching out the BED file we can then generate a whole genome high
coverage read file.

### Behavior at region edges

If the start of a region (defined in reference coordinates) falls in the middle of a deletion, the actual start of
the region is shifted to be the first base before the deletion. This allows reads generated to be consistent with
a VCF that has been cut using that BED file. If the end of a region cuts a deletion the actual end
is shifted to the first base after the deletion so reads capture that deletion. **If there is a SNP on
that base, it will be included in the reads and this is not consistent with a VCF cut with that BED file** It is best
to ensure this second case does not happen by altering the BED file to either stop before the deletion, or include
the neighboring SNP.

Generation of Pos and CIGAR strings
-----------------------------------

(The implementation of this algorithm is in `mitty.simulation.rpc`)

Variants are applied to the reference sequence to generate a sequence of nodes.
Each node carries a sequence fragment and we can concatenate the sequence fragments
from neighboring nodes to create a read. More importantly, each node carries a set
of metadata that allows us to

   - Indicate the correct POS and CIGAR value for that read
   - Indicate the variants carried by the read

Example
=======

Reference sequence

```
ATGACGTATCCAAGGAGGCGTTACC
1234567890123456789012345
```

Variants

```
##fileformat=VCFv4.1
##contig=<ID=tiny,length=51>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	g0_s0
1	5	.	C	T	100	PASS	.	GT	0|1
1	8	.	A	ATTT	100	PASS	.	GT	0|1
1	11	.	CAA	C	100	PASS	.	GT	0|1
1	14	.	G	T	100	PASS	.	GT	1|0
1	20	.	GTT	G	100	PASS	.	GT	1|1
```

For chrom copy 0|1

```
+---------------+---------+-----------------------------------------------------+-------+
|     REF       |  ALT    |                      Node data                              |
+-------+-------+  0|1    +------------+-----------+----------+--------+--------+-------+
|  POS  |  SEQ  |         |    seq     |    pr     | CIGAR OP | OP LEN |   ps   |v [1]  |
+=======+=======+=========+============+===========+==========+========+========+=======+
|   1   |   A   |    A    |    ATGA    |    1      |    M     |    4   |   1    |       |
+-------+-------+---------+            |           |          |        |        |       |
|   2   |   T   |    T    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|   3   |   G   |    G    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|   4   |   A   |    A    |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|   5   |   C   |    T    |     T      |    5      |    X     |    1   |   5    |   0   |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|   6   |   G   |    G    |    GTA     |    6      |    M     |    3   |   6    |       |
+-------+-------+---------+            |           |          |        |        |       |
|   7   |   T   |    T    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|   8   |   A   |    A    |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|               |    T    |    TTT     |    9      |    I     |    3   |   9    |   3   |
|               +---------+            |           |          |        |        |       |
|               |    T    |            |           |          |        |        |       |
|               +---------+            |           |          |        |        |       |
|               |    T    |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|   9   |   T   |    T    |    TCC     |    9      |    M     |    3   |   12   |       |
+-------+-------+---------+            |           |          |        |        |       |
|  10   |   C   |    C    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  11   |   C   |    C    |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|  12   |   A   |         |            |   14 [2]  |    D     |    2   | 14 [3] |  -2   |
+-------+-------+         |            |           |          |        |        |       |
|  13   |   A   |         |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|  14   |   G   |    G    |  GGAGGCG   |    14     |    M     |    7   |   15   |       |
+-------+-------+---------+            |           |          |        |        |       |
|  15   |   G   |    G    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  16   |   A   |    A    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  17   |   G   |    G    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  18   |   G   |    G    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  19   |   C   |    C    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  20   |   G   |    G    |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|  21   |   T   |         |            |    23     |    D     |   2    |   21   |  -2   |
+-------+-------+         |            |           |          |        |        |       |
|  22   |   T   |         |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
|  23   |   A   |    A    |    ACC     |    23     |    M     |   3    |   22   |       |
+-------+-------+---------+            |           |          |        |        |       |
|  24   |   C   |    C    |            |           |          |        |        |       |
+-------+-------+---------+            |           |          |        |        |       |
|  25   |   C   |    C    |            |           |          |        |        |       |
+-------+-------+---------+------------+-----------+----------+--------+--------+-------+
```

[1] Size of variant
[2] Not needed during read generation computations
[3] Need to set this to the ps of the last base before the deletion to make the CIGAR string generation work


Reference sequence::

  ATGACGTATCCAAGGAGGCGTTACC
  1234567890123456789012345

Variants::

  ##fileformat=VCFv4.1
  ##contig=<ID=tiny,length=51>
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	g0_s0
  1	5	.	C	T	100	PASS	.	GT	0|1
  1	8	.	A	ATTT	100	PASS	.	GT	0|1
  1	11	.	CAA	C	100	PASS	.	GT	0|1
  1	14	.	G	T	100	PASS	.	GT	1|0
  1	20	.	GTT	G	100	PASS	.	GT	1|1


For chrom copy 1|0

+---------------+---------+-----------------------------------------------------------------+
|     REF       |  ALT    |                      Node data                                  |
+-------+-------+  1|0    +----------------+-----------+----------+--------+--------+-------+
|  POS  |  SEQ  |         |    seq         |    pr     | CIGAR OP | OP LEN |   ps   |   v   |
+=======+=======+=========+================+===========+==========+========+========+=======+
|   1   |   A   |    A    | ATGACGTATCCAA  |    1      |    M     |   13   |   1    |       |
+-------+-------+---------+                |           |          |        |        |       |
|   2   |   T   |    T    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   3   |   G   |    G    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   4   |   A   |    A    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   5   |   C   |    C    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   6   |   G   |    G    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   7   |   T   |    T    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   8   |   A   |    A    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|   9   |   T   |    T    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  10   |   C   |    C    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  11   |   C   |    C    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  12   |   A   |    A    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  13   |   A   |    A    |                |           |          |        |        |       |
+-------+-------+---------+----------------+-----------+----------+--------+--------+-------+
|  14   |   G   |    T    |       T        |    14     |    X     |   1    |   14   |   0   |
+-------+-------+---------+----------------+-----------+----------+--------+--------+-------+
|  15   |   G   |    G    |     GAGGCG     |    15     |    M     |   6    |   15   |       |
+-------+-------+---------+                |           |          |        |        |       |
|  16   |   A   |    A    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  17   |   G   |    G    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  18   |   G   |    G    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  19   |   C   |    C    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  20   |   G   |    G    |                |           |          |        |        |       |
+-------+-------+---------+----------------+-----------+----------+--------+--------+-------+
|  21   |   T   |         |                |    23     |    D     |   2    |   20   |   -2  |
+-------+-------+         |                |           |          |        |        |       |
|  22   |   T   |         |                |           |          |        |        |       |
+-------+-------+---------+----------------+-----------+----------+--------+--------+-------+
|  23   |   A   |    A    |      ACC       |    23     |    M     |   3    |   21   |       |
+-------+-------+---------+                |           |          |        |        |       |
|  24   |   C   |    C    |                |           |          |        |        |       |
+-------+-------+---------+                |           |          |        |        |       |
|  25   |   C   |    C    |                |           |          |        |        |       |
+-------+-------+---------+----------------+-----------+----------+--------+--------+-------+


Creating the node sequence
==========================
Say that the reference sequence is ``rseq``, the current cursor on the reference sequence is at
``P_ref`` and the current variant is ``(P_v, ref, alt)``. We will create two nodes to
add to the sequence according to the following rules.

+------------+----------------------------------------------------------------------------+
|            |                               Node data                                    |
+------------+---------------------------+-----------+---------+---------+----------------+
| Variant    |  seq                      |           |         |  CIGAR  |   op           |
| type       |                           |   pr      |   ps    |  op     |   len          |
+============+===========================+===========+=========+=========+================+
|            | ``rseq[P_ref .. P_v - 1]``| ``P_ref`` |  ``P*`` |   M     |  len(seq)      |
| SNP        +---------------------------+-----------+---------+---------+----------------+
|            | ``alt``                   | ``P_v``   |  ``P*`` |   X     |  len(seq)      |
+------------+---------------------------+-----------+---------+---------+----------------+
|            | ``rseq[P_ref .. P_v]``    | ``P_ref`` |  ``P*`` |   M     |  len(seq)      |
| INS        +---------------------------+-----------+---------+---------+----------------+
|            | ``alt[2 ... ]``           | ``-``     |  ``P*`` |   I     |  len(seq)      |
+------------+---------------------------+-----------+---------+---------+----------------+
|            | ``rseq[P_ref .. P_v]``    | ``P_ref`` |  ``P*`` |   M     |  len(seq)      |
| DEL        +---------------------------+-----------+---------+---------+----------------+
|            |                           | ``-``     |  ``-``  |   D     |  len(alt) - 1  |
+------------+---------------------------+-----------+---------+---------+----------------+

Here ``p*`` is the cumulative count of sample sequence bases at the start of the node. It
is the cumulative sum of len(seq) for all nodes except D nodes

If an M-node (the first node) has no sequence (i.e. op len = 0) it is omitted

Generating reads and alignment metadata
=======================================

A read starting at position ``P_samp`` with length ``r_len`` is generated as follows:

- The first node with ``ps`` less than ``P_samp`` is the starting node: ``N[1]``
- The node with ``ps + len(seq)`` less than ``P_samp + r_len`` is the end node: ``N[m]``
- By this definition, a D node can't be a start or end node.
- We concatenate the node sequence data together as follows: for ``N[2] ... N[m-1]`` we
  concatenate all the sequence information, and for the start and end ndoes N[1] and N[m]
  we slice the sequence depending on where ``P_samp`` starts and ``P-samp + r_len`` ends
- We similarly concatenate the CIGAR ops, prefixing them with the ``op len`` value,
  similarly curtailing for ``N[1]`` and ``N[m]``
- We also create a list of variants the read carries by concatenating all the non-M nodes
  with the ``op_len`` as variant length and sign being determined by type of variant (D = -1,
  I = +1). SNPs are an exception as their variant length is set as 0.