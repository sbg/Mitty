Read Generation Algorithms
++++++++++++++++++++++++++


Reducing bias in read location
------------------------------

- Read positions and lengths over the whole genome are pre-generated and stored in an array

  - This can be done repeatedly if main memory is not sufficient to store all the required data
- Chromosomes are broken up into small segments (lowerbound on size is governed by the maximum length of read)
- These segments are fed in a scrambled order to parallel threads

  - the read list is checked to find appropriate reads that can be generated for the segment
  - the templates/reads are fed to the output in the order they return.
  - The multi-threading creates some uncontrolled randomness to the order and the scrambling of the segments
    adds some controlled randomness to the order of the reads generated


Generation of Pos and CIGAR strings
-----------------------------------

::

        12345678901
  ref - ATGACGTATCC

      POS  REF  ALT
      3    G    T
      5    C    CTT
      7    TAT  T


              12  3  45      67  89  01
       ref -  AT  G  AC      GT  AT  CC

      samp -  AT  T  AC  TT  GT  |   CC
        pr -  12  3  45      67      01
        ps -  12  3  45  67  89      01
  cigar op -  M   X  M   I   M   D   M
    op len -  2   1  2   2   2   2   2

              ##--#--##--##           -->  POS: 1, CIGAR: 2=1X2=2I
               #--#--##--##--#         -->  POS: 2, CIGAR: 1=1X2=2I1M
                  #--##--##--##         -->  POS: 3, CIGAR: 1X2=2I2M
                     ##--##--##------#   -->  POS: 4, CIGAR: 2=2I2M2D1M
                      #--##--##------##   -->  POS: 5, CIGAR: 1=2I2M2D2M
                         ##--##------##    ---> POS: 6, CIGAR: 2I2M2D2M


Variants are applied to the reference sequence to generate a sequence of nodes.
Each node carries a sequence fragment and we can concatenate the sequence fragments
from neighboring nodes to create a read. More importantly, each node carries a set
of metadata that allows us to

   - Indicate the correct POS and CIGAR value for that read
   - Indicate the variants carried by the read


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

If an M-node (the first node) as no sequence (i.e. op len = 0) it is omitted

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


