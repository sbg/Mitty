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


class Node(object):
  __slots__ = ('ps', 'pr', 'cigarop', 'oplen', 'seq', 'v')

  def __init__(self, ps, pr, cigarop, oplen, seq, v):
    self.ps = ps
    self.pr = pr
    self.cigarop = cigarop
    self.oplen = oplen
    self.seq = seq
    self.v = v  # The variant size code. None for reference matching nodes

  def tuple(self):
    return self.ps, self.pr, self.cigarop, self.oplen, self.seq, self.v

  def __repr__(self):
    return self.tuple().__repr__()

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.tuple() == other.tuple()
    else:
      return self.tuple() == other

  def __ne__(self, other):
    return not self.__eq__(other)


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
  for v in vcf.itertuples():
    if v.GT & chrom_copy:  # This variant exists on this chrom copy
      if v.POS < ref_pos: continue  # We are starting from a later position
      new_nodes, samp_pos, ref_pos = create_nodes(ref_seq, samp_pos, ref_pos, v)
      nodes += new_nodes

  delta = ref_pos - ref_start_pos
  if delta <= len(ref_seq):  # Last part of sequence, needs an M
    nodes.append(Node(samp_pos, ref_pos, '=', len(ref_seq) - delta, ref_seq[delta:], None))

  return nodes


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
      raise RuntimeError("Complex variants present in VCF. Please filter or refactor these.")
    if l_a > 1:  # INS
      return 'I', l_a - l_r
    else:  # SNP
      return 'X', 0
  else:
    if l_a > 1:
      raise RuntimeError("Complex variants present in VCF. Please filter or refactor these.")
    return 'D', l_r - l_a


def snp(ref_seq, samp_pos, ref_pos, v, vl):
  nodes = []
  delta = v.POS - ref_pos
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append(Node(samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - 1:v.POS - 1], None))  # M node
    ref_pos = v.POS
    samp_pos += delta
  #               ps         pr    op oplen seq
  nodes.append(Node(samp_pos, ref_pos, 'X', 1, v.ALT, 0))
  ref_pos += 1
  samp_pos += 1
  return nodes, samp_pos, ref_pos


def insertion(ref_seq, samp_pos, ref_pos, v, vl):
  nodes = []
  delta = v.POS + 1 - ref_pos  # This is INS, so we include the first base in the M segment
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append(Node(samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - 1:v.POS], None))  # M node
    samp_pos += delta

  ref_pos = v.POS + 1  # The next ref pos is the ref base just after the insertion
  #                ps        pr     op  oplen  seq
  nodes.append(Node(samp_pos, ref_pos, 'I', vl, v.ALT[1:], vl))
  samp_pos += vl
  return nodes, samp_pos, ref_pos


def deletion(ref_seq, samp_pos, ref_pos, v, vl):
  nodes = []
  delta = v.POS + 1 - ref_pos  # This is DEL, so we include the first base in the M segment
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append(Node(samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - 1:v.POS], None))  # M node
    samp_pos += delta

  ref_pos = v.POS + 1 + vl  # The next ref pos is the ref base just after the deletion
  #                    ps           pr     op  oplen    seq
  nodes.append(Node(samp_pos - 1, ref_pos, 'D', vl, '', -vl))
  return nodes, samp_pos, ref_pos


def get_begin_end_nodes(pl, ll, nodes):
  """Given a list of read positions and lengths return us a list of start and end nodes

  :param pl:  should be np.array so we can sum pl and ll
  :param ll:             "
  :param nodes:
  :return: n0, n1 each is a list of nodes
  """
  ps = np.array([n.ps if n.cigarop != 'D' else n.ps + 1 for n in nodes], dtype=np.uint64)
  # D nodes .ps values are set to last base before deletion. We increment this by one so
  # we can get proper start/stop node computation
  return ps.searchsorted(pl, 'right') - 1, ps.searchsorted(pl + ll - 1, 'right') - 1


def generate_read(p, l, n0, n1, nodes):
  """

  :param p: Start position of read in sample coordinates
  :param l: Length of read
  :param n0: starting node
  :param n1: ending node
  :param nodes: as returned by create_node_list
  :return: (pos, cigar, v_list, seq)
           v_list = [-d, +i, 0] -> a list of ints indicating size of variants carried by read
  """
  v_list = [n.v for n in nodes[n0:n1 + 1] if n.v is not None]
  cigar = [str(min(p + l - n.ps, n.oplen) - max(0, p - n.ps)) + n.cigarop for n in nodes[n0:n1 + 1]]
  seq = [n.seq[max(0, p - n.ps):min(p + l - n.ps, n.oplen)] for n in nodes[n0:n1 + 1]]

  if nodes[n0].cigarop == 'I':
    if n0 == n1:
      # Special case - read is from inside a long insertion
      # We want to pile up reads from a long insertion at the start of the insertion
      # We need to override the CIGAR too
      pos = nodes[n0].pr - 1
      cigar = ['>{}:{}I'.format(p - nodes[n0].ps, l)]
    else:
      pos = nodes[n0].pr
  else:
    pos = p - nodes[n0].ps + nodes[n0].pr

  return pos, ''.join(cigar), v_list, ''.join(seq)