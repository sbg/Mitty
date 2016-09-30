"""
Algorithms:

Read location selection:

Given p and N generate a coverage of pN = C by taking a section of the
reference sequence + VCF and running the read generation algorithm N times with
a probability of success p each time.

Pos and CIGAR computation:

The sample sequence is represented as a list of nodes

(vertex, (vertex, vertex))
   |         |       |
   |         |       --------- vertex for chrom copy 2
   |         ----------------- vertex for chrom copy 1
   --------------------------- current vertex

vertex = (pos, cigarop, seq)
           |       |     |
           |       |     ----- sequence. None for deletions
           |       ----------- cigarop  'M', 'I' 'D' or 'X'
           ------------------- pos on ref branch for start of node
                               None for D and I

offset = (vertex, offset)
             |      |
             |      ---------- relative to start of vertex.
             ----------------- can't be a D vertex




"""

__qname_format__ = '@read_serial|chrom|copy|strand|pos|rlen|cigar|vs1,vs2,...|strand|pos|rlen|cigar|vs1,vs2,...'
__qname_format_details__ = """
@read_serial|chrom|copy|strand|pos|rlen|cigar|vs1,vs2,...|strand|pos|rlen|cigar|vs1,vs2,...
    |          |     |    |     |    |    |        |         |                      |
 unique        |     |    |     | read    |        |         ---- repeated for ------
 number for    |     |    |     | len     |        |          other read in template
 template      |     |    |     |         |        |
               |     |    |     |     cigar    comma separated
      chrom read     |    |     |              list of sizes of
  was taken from     |    |     |              variants this read
       One based     |    |     |              covers
                     |    |     |
                     |    |     |
    copy of chrom read    |     |
  was taken from (1,2)    |     |
                          |     |
         forward strand (0)     |
      or reverse strand (1)     |
                                |
                      pos of read
                        One based

The chrom and pos are one based to make comparing qname info in genome browser easier

For reads from inside a long insertion the CIGAR has the following format:

  '>p:nI'

where:

 '>' is the unique key that indicates a read inside a long insertion
 'p' is how many bases into the insertion branch the read starts
 'n' is simply the length of the read
"""
import numpy as np


def create_node_list(ref_seq, ref_start_pos, chrom_copy, vcf):
  """Apply the apropriate variants (i.e. those that are on this copy of the
  chromosome) to the reference sequence and create the sequence of nodes which
  other functions use to generate reads and alignment metadata.

  :param ref_seq: Reference sequence
  :param ref_start_pos: position where this part of the reference sequence starts from. 1 indexed
  :param chrom_copy: b01, b10, b11
  :param vcf: variant information, same format as VCF file
              Should be from correct chrom.
              Restricting the list of variants to just those
              covering the passed sequence speeds up processing

  :return:

  """
  samp_pos, ref_pos = 1, ref_start_pos  # 1 indexed
  # These are both 1 indexed.
  # ref_pos is relative to whole sequence, samp_pos is relative to this fragment of the expanded sequence

  nodes = []
  for _, v in vcf.itertuples():
    if v.GT & chrom_copy:  # This variant exists on this chrom copy
      new_nodes, samp_pos, ref_pos = create_nodes(ref_seq, samp_pos, ref_pos, v)
      nodes += new_nodes

  delta = ref_pos - ref_start_pos
  if delta <= len(ref_seq):  # Last part of sequence, needs an M
    nodes.append((samp_pos, ref_pos, 'M', len(ref_seq) - delta, ref_seq[delta:]))

  return np.array(nodes,
                  dtype=[('ps', np.uint32),
                         ('pr', np.uint32),
                         ('cigarop', np.ubyte),
                         ('oplen', np.uint32),
                         ('seq', str)])


def create_nodes(ref_seq, samp_pos, ref_pos, v):
  vt, vl = get_v_type(v)
  if vt == 'X':
    return snp(ref_seq, samp_pos, ref_pos, v, vl)
  elif vt == 'I':
    return insertion(ref_seq, samp_pos, ref_pos, v, vl)
  else:
    return deletion(ref_seq, samp_pos, ref_pos, v, vl)


def get_v_type(v):
  """Return us the variant type and size

  :param v:
  :return:
  """
  l_r, l_a = len(v.REF), len(v.ALT)
  # Figure out what kind of a variant it is
  if l_r <= l_a:  # SNP or INS
    if l_r > 1:
      raise RuntimeError('Complex variants present in VCF. Please filter or refactor these.')
    if l_a > 1:  # INS
      return 'I', l_a - l_r
    else:  # SNP
      return 'X', 0
  else:
    if l_a > 1:
      raise RuntimeError('Complex variants present in VCF. Please filter or refactor these.')
    return 'D', l_a - l_r  # Will be negative


def snp(ref_seq, samp_pos, ref_pos, v, vl):
  nodes = []
  delta = v.POS - ref_pos
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append((samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - 1:v.POS - 1]))  # M node
    ref_pos = v.POS
    samp_pos += delta
  #               ps         pr     op  oplen  seq
  nodes.append((samp_pos, ref_pos, 'X', 0, v.ALT))
  ref_pos += 1
  samp_pos += 1
  return nodes, samp_pos, ref_pos


def insertion(ref_seq, samp_pos, ref_pos, v, vl):
  nodes = []
  delta = v.POS + 1 - ref_pos  # This is INS, so we include the first base in the M segment
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append((samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - 1:v.POS]))  # M node
    ref_pos = v.POS + 1  # The next ref pos is the ref base just after the insertion
    samp_pos += delta

  #                ps        pr     op  oplen  seq
  nodes.append((samp_pos, ref_pos, 'I', vl, v.ALT[1:]))
  samp_pos += vl
  return nodes, samp_pos, ref_pos


def deletion(ref_seq, samp_pos, ref_pos, v, vl):
  nodes = []
  delta = v.POS + 1 - ref_pos  # This is DEL, so we include the first base in the M segment
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append((samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - 1:v.POS]))  # M node
    ref_pos = v.POS + 1  # The next ref pos is the ref base just after the insertion
    samp_pos += delta

  #                ps        pr     op  oplen  seq
  nodes.append((samp_pos, ref_pos, 'I', vl, v.ALT[1:]))
  samp_pos += vl
  return nodes, samp_pos, ref_pos


def small_section(ref, variants, p, N):
  pass