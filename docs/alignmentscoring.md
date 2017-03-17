Alignment scoring
=================
The read is scored using the metric `d_err` which states a difference in base pairs between the aligned position of the read and the correct position. Positive `d_err` values mean the read is aligned too far to the right.

Under the "strict" scoring scheme the read POS has to exactly match the correct POS. When a linient scoring scheme is used, if the aligner has soft or hardclipped any bases at the start of the read, the comparison point is moved to the right by that many bases, effectively checking if the aligner has correctly placed the read after the softclip. 

By running a contrast analyses with these two settings it is possible to uncover problems with the soft-clipping algorithms used by aligners.

Reads from completely inside insertions, which are normally left unmapped by a linear aligner, can be mapped by a graph aligner if the insertion is known. If the aligner piles up these reads at the start of the insertion point, an reassembly variant caller like Haplotype caller can retrieve the insertion. For this reason, such reads are scored relative to the insertion location.