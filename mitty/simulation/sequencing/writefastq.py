"""Utility functions to write simulated reads into a FASTQ file. This takes care of writing
the sidecar file for all templates whose qname length exceeds the 254 character limit set
in the sam spec."""
import time
import logging
from collections import namedtuple

from numpy import base_repr

logger = logging.getLogger(__name__)
__process_stop_code__ = 'SETECASTRONOMY'

__max_qname_len__ = 240
__qname_format_details__ = """
      @index|sn|chrom:copy|strand:pos:cigar:v1,v2,...:MD|strand:pos:cigar:v1,v2,...:MD*

index:  unique code for each template. If you must know, the index is simply a monotonic counter
        in base-36 notation
sn:     sample name. Useful if simulation is an ad-mixture of sample + contaminants or multiple samples
chrom:  chromosome id of chromosome the read is taken from
copy:   copy of chromosome the read is taken from (0, 1, 2, ...) depends on ploidy
strand: 0: forward strand, 1: reverse strand
pos:    Position of first (left-most) reference matching base in read (Same definition as for POS in BAM)
        One based
        For reads coming from completely inside an insertion this is, however, the POS for the insertion
cigar:  CIGAR string for correct alignment.
        For reads coming from completely inside an insertion, however, the CIGAR string is:
        '>p+nI' where:
           '>' is the unique key that indicates a read inside a long insertion
           'p' is how many bases into the insertion branch the read starts
           'n' is simply the length of the read
v1,v2,..: Comma separated list of variant sizes covered by this read
            0 = SNP
            + = INS
            - = DEL
MD:     Read corruption MD tag. This is empty for perfect reads, but is filled with an MD formatted
        string for corrupted reads. The MD string is referenced to the original, perfect read, not
        a reference sequence.

Notes:

- The alignment information is repeated for every read in the template.
  The example shows what a PE template would look like.
- The pos value is are one based to make comparing qname info in genome browser easier
- qnames longer than N characters are stored in a side-car file alongside the simulated FASTQs.
  The qname in the FASTQ file itself is truncated to N characters.

  Nominally, N=254 according to the SAM spec, but due to bugs in some versions of htslib
  this has been set shorter to {}.
  A truncated qname is detected when the last character is not *
""".format(__max_qname_len__)


def writer(fastq1_out, side_car_out, fastq2_out=None, data_queue=None, max_qname_len=__max_qname_len__):
  """Write templates to file

  :param fastq1_out: Name of FASTQ1
  :param side_car_out: Name of side car file for long qnames
  :param fastq2_out: If paired end, name of FASTQ2
  :param data_queue: multiprocessing queue
  :param max_qname_len: Send qnames longer than this to the overflow file

    The data format is as follows:
      (
        idx,  - if this is None then we create an index afresh
        sample_name,
        chrom,
        copy,
        (
          (strand, pos, cigar, (v1,v2,...), MD, seq, qual)
        ...  [repeated as for as many reads in this template]
        )
      )

  :return:
  """
  t0 = time.time()

  cnt = -1
  fastq_l, side_car_fp = [open(fastq1_out, 'w')], open(side_car_out, 'w')
  if fastq2_out is not None: fastq_l += [open(fastq2_out, 'w')]
  for cnt, template in enumerate(iter(data_queue.get, __process_stop_code__)):
    # @index|sn|chrom:copy|
    qname = '@{}|{}|{}:{}|'.format(template[0] or base_repr(cnt, 36), *template[1:4])
    # strand:pos:cigar:v1,v2,...:MD|strand:pos:cigar:v1,v2,...:MD*
    qname += '|'.join('{}:{}:{}:{}:{}'.format(*r[:3], str(r[3])[1:-1].replace(' ', ''), r[4]) for r in template[4]) + '*'

    if len(qname) > max_qname_len:
      side_car_fp.write(qname + '\n')
      qname = qname[:max_qname_len]

    for fp, r in zip(fastq_l, template[4]):
      fp.write('{}\n{}\n+\n{}\n'.format(qname, r[5], r[6]))

  for fp in fastq_l:
    fp.close()

  t1 = time.time()
  logger.debug('Writer finished: {} templates in {:0.2f}s ({:0.2f} t/s)'.format(cnt + 1, t1 - t0, (cnt + 1) / (t1 - t0)))


ri = namedtuple('ReadInfo',
                ['index', 'sample', 'chrom', 'cpy', 'strand', 'pos', 'cigar', 'special_cigar', 'v_list', 'md'])


def parse_qname(qname, long_qname_table=None):
  """Given a Mitty qname return us the POS and CIGAR as we would put in a BAM. There is also a special_cigar
  which is set for reads completely inside long insertions

  @index|sn|chrom:copy|strand:pos:cigar:v1,v2,...:MD|strand:pos:cigar:v1,v2,...:MD|

  :param qname:
  :param long_qname_table: If present, is a map of qname index and qname for just the long qnames
  :return: pos, cigar, special_cigar
  """
  def _parse_(_cigar, _v_list):
    """
    Parse cigar to extract special_cigar if needed
    Parse v_list

    :param _cigar:
    :param _v_list:
    :return:
    """
    if _cigar[0] == '>':  # This is a special_cigar for a read from inside an insertion
      _special_cigar = _cigar
      _cigar = _cigar.split('+')[-1]
    else:
      _special_cigar = None

    return _cigar, _special_cigar, [int(v) for v in _v_list.split(',') if v is not '']

  def _split_(_r):
    strand, pos, cigar, v_list, md = _r.split(':')
    return (int(strand), int(pos)) + _parse_(cigar, v_list) + (md,)

  if qname[-1] != '*':  # Truncated qname
    if long_qname_table is None:
      raise ValueError('Long qname with no table lookup')  # It's the caller's responsibility to handle this error
    _qname = long_qname_table.get(qname.split('|', 1)[0], None)
    if _qname is None:
      raise ValueError('Long qname with no table lookup: {}, {}'.format(qname.split('|', 1)[0], qname))
      # It's the caller's responsibility to handle this error
    qname = _qname

  d = qname[:-1].split('|')  # Strip out the terminal '*'
  serial, sample = d[:2]
  chrom, cpy = d[2].split(':')
  cpy = int(cpy)
  return [
    ri(serial, sample, chrom, cpy, *_split_(r))
    for r in d[3:]
  ]


# This assumes that the longer qnames are relatively rare, so this sidecar file is relatively small
# and system memory is sufficient to store all such exceptions at once
def load_qname_sidecar(sidecar_fn):
  def parse_sidecar_file(_fn):
    with open(_fn, 'r') as fp:
      for ln in fp:
        yield ln.split('|', 1)[0][1:], ln[1:-1]
  return {
    k: v
    for k, v in parse_sidecar_file(sidecar_fn)
  }