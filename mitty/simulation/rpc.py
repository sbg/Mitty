"""Algorithms for read generation"""
import numpy as np


class Node(object):
  __slots__ = ('ps', 'pr', 'cigarop', 'oplen', 'seq', 'v')

  def __init__(self, ps, pr, cigarop, oplen, seq):
    self.ps = ps
    self.pr = pr
    self.cigarop = cigarop
    self.oplen = oplen
    self.seq = seq
    # The variant size code. None for reference matching nodes
    self.v = {
      '=': None,
      'X': 0,
      'I': oplen,
      'D': -oplen
    }[cigarop]

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


def create_node_list(ref_seq, ref_start_pos, vl):
  """Apply the appropriate variants (i.e. those that are on this copy of the
  chromosome) to the reference sequence and create the sequence of nodes which
  other functions use to generate reads and alignment metadata.

  :param ref_seq: Reference sequence
  :param ref_start_pos: position where this part of the reference sequence starts from. 1 indexed
  :param vl: list of Variant objects in order
  :return:

  """
  samp_pos, ref_pos = ref_start_pos, ref_start_pos  # 1 indexed
  # These are both 1 indexed.
  # ref_pos is relative to whole sequence, samp_pos is relative to this fragment of the expanded sequence

  nodes = []
  for v in vl:
    if v.pos < ref_pos: continue  # We are starting from a later position
    new_nodes, samp_pos, ref_pos = create_nodes(ref_seq, samp_pos, ref_pos, v, ref_start_pos)
    nodes += new_nodes

  offset = ref_pos - ref_start_pos
  if offset <= len(ref_seq):  # Last part of sequence, needs an M
    nodes.append(Node(samp_pos, ref_pos, '=', len(ref_seq) - offset, ref_seq[offset:]))

  return nodes


def create_nodes(ref_seq, samp_pos, ref_pos, v, ref_start_pos):
  if v.cigarop == 'X':
    return snp(ref_seq, samp_pos, ref_pos, v, ref_start_pos)
  elif v.cigarop == 'I':
    return insertion(ref_seq, samp_pos, ref_pos, v, ref_start_pos)
  else:
    return deletion(ref_seq, samp_pos, ref_pos, v, ref_start_pos)


def snp(ref_seq, samp_pos, ref_pos, v, ref_start_pos):
  nodes = []
  delta = v.pos - ref_pos
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append(Node(samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - ref_start_pos:v.pos - ref_start_pos]))  # M node
    ref_pos = v.pos
    samp_pos += delta
  #               ps         pr    op oplen seq
  nodes.append(Node(samp_pos, ref_pos, 'X', 1, v.alt))
  ref_pos += 1
  samp_pos += 1
  return nodes, samp_pos, ref_pos


def insertion(ref_seq, samp_pos, ref_pos, v, ref_start_pos):
  nodes = []
  delta = v.pos + 1 - ref_pos  # This is INS, so we include the first base in the M segment
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append(Node(samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - ref_start_pos:v.pos + 1 - ref_start_pos]))  # M node
    samp_pos += delta

  ref_pos = v.pos + 1  # The next ref pos is the ref base just after the insertion
  #                ps        pr     op  oplen  seq
  nodes.append(Node(samp_pos, ref_pos, 'I', v.oplen, v.alt[1:]))
  samp_pos += v.oplen
  return nodes, samp_pos, ref_pos


def deletion(ref_seq, samp_pos, ref_pos, v, ref_start_pos):
  nodes = []
  delta = v.pos + 1 - ref_pos  # This is DEL, so we include the first base in the M segment
  if delta > 0:  # Need to make an M node
    #               ps         pr     op   oplen          seq
    nodes.append(Node(samp_pos, ref_pos, '=', delta, ref_seq[ref_pos - ref_start_pos:v.pos + 1 - ref_start_pos]))  # M node
    samp_pos += delta

  ref_pos = v.pos + 1 + v.oplen  # The next ref pos is the ref base just after the deletion
  #                    ps           pr     op  oplen    seq
  nodes.append(Node(samp_pos - 1, ref_pos, 'D', v.oplen, ''))
  return nodes, samp_pos, ref_pos


def get_begin_end_nodes(pl, ll, nodes):
  """Given a list of read positions and lengths return us a list of start and end nodes

  :param pl:  Positions of reads. should be np.array so we can sum pl and ll
  :param ll:  Lengths of reads.           "
  :param nodes:
  :return: nse 2 x N np.array (N = len(pl)
  """
  ps = np.array([n.ps if n.cigarop != 'D' else n.ps + 1 for n in nodes], dtype=np.uint64)
  # D nodes .ps values are set to last base before deletion. We increment this by one so
  # we can get proper start/stop node computation
  return [ps.searchsorted(pl, 'right') - 1, ps.searchsorted(pl + ll - 1, 'right') - 1]


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
  pos = nodes[n0].pr
  cigar = [str((min(p + l - n.ps, n.oplen) - max(0, p - n.ps)) if n.cigarop != 'D' else n.oplen) + n.cigarop for n in nodes[n0:n1 + 1]]
  v_list = [n.v for n in nodes[n0:n1 + 1] if n.v is not None]
  seq = [n.seq[max(0, p - n.ps):min(p + l - n.ps, n.oplen)] for n in nodes[n0:n1 + 1]]

  if nodes[n0].cigarop == 'I':
    if n0 == n1:
      # Special case - read is from inside a long insertion
      # We want to pile up reads from a long insertion at the start of the insertion, which happens automatically
      # Now we need to override the CIGAR
      cigar = ['>{}+{}I'.format(p - nodes[n0].ps, l)]
  else:
    pos = p - nodes[n0].ps + nodes[n0].pr

  return pos, ''.join(cigar), v_list, ''.join(seq)