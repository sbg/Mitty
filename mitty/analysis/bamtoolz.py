"""
Takes care of reading from BAM files and scattering (parallelizing) pipelines

See docs/bamtoolz.md for algorithm details"""
from collections import OrderedDict
import time
from multiprocessing import Process, Queue
import queue

import pysam
import cytoolz.curried as cyt

import logging
logger = logging.getLogger(__name__)


def get_seq_dict(bam_fname):
  """Return us a dictionary converting reference_id to sequence name

  :param bam_fname:
  :return:
  """
  fp = pysam.AlignmentFile(bam_fname)
  contig_dict = {n: cn['SN'] for n, cn in enumerate(fp.header['SQ'])}
  contig_dict[-1] = '*'
  return contig_dict


def read_bam_st(bam_fname):
  """Vanilla single thread single read iterator

  :param bam_fname:
  :return: iterator over single read tuples (read,)
  """
  for read in pysam.AlignmentFile(bam_fname).fetch(until_eof=True):
    if read.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    yield (read,)


def read_bam_paired_st(bam_fname):
  """Vanilla single thread paired read iterator

  :param bam_fname:
  :return: iterator over paired read tuples (read1, read2)
  """
  singles = {}
  for read in pysam.AlignmentFile(bam_fname).fetch(until_eof=True):
    if read.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
    key = read.qname[:20]  # Is this enough?
    if key not in singles:
      singles[key] = read
    else:  # You complete me
      yield (read, singles[key]) if read.is_read1 else (singles[key], read)
      del singles[key]
  if len(singles):
    logger.error('{} unpaired reads left over!'.format(len(singles)))
    logger.error(singles.keys())


def read_iter(fp, contig_q):
  """Returns read objects from contigs until someone passes None as a contig

  :param fp: BAM file pointer (pysam.AlignmentFile)
  :param contig_q: a queue into which we put contig information
                   (contig, eof_true) - eof_true is set if this is
                   the last non-empty contig and we want to pull out
                   all the trailing unmapped reads right after the contig
  :return: a generator
  """
  for contig in iter(contig_q.get, None):
    logger.debug(contig[0])
    for read in fp.fetch(contig[0]):
      if read.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
      yield read
    if contig[1]:  # Now want the trailing reads - fp is positioned just before them
      for read in fp.fetch(until_eof=contig[1]):
        if read.flag & 0b100100000000: continue  # Skip supplementary or secondary alignments
        yield read


def unpaired_read_iter(fp, contig_q):
  """Almost identical to read_iter, except it returns tuples (read,)
  This enables us to write processing code that operates both on single
  reads as well as pairs since they both come as tuples

  :param fp: BAM file pointer (pysam.AlignmentFile)
  :param contig_q: a queue into which we put contig information
                   (contig, eof_true) - eof_true is set if this is
                   the last non-empty contig and we want to pull out
                   all the trailing unmapped reads
  :return: a generator that yields (read,) tuples
  """
  for read in read_iter(fp, contig_q):
    yield (read,)


def paired_read_iter(fp, contig_q, singles_q, max_singles=1000,
                     is_singles_mixer=False, single_src_cnt=None):
  """

  :param fp:                pysam.AlignmentFile()

  :param contig_q:          Messages are of the form (ref, eof)
                            ref is the name of the contig to fetch
                            eof is T/F and indicates whether we should fetch till eof
                            the sender should set this to T only if this is the last
                            non-empty contig and is followed by the unmapped reads
                            that sit at the end of the file.

  :param singles_q:         Messages are SAM strings of reads converted using tostring()
                            reads are recieved and mixed in with the from-disk stream
                            if this is the singles mixer, else, unpaired reads are
                            sent to this Q

  :param max_singles:       When we have these many singles, start passing then to the
                            singles mixer

  :param is_singles_mixer:  Set True if this is also the "singles mixer" that
                            receives unpaired reads from other workers

  :param single_src_cnt:    How many processes are out there sending singles?
                            Used if this is a singles mixer

  :return: a generator that yields paired read tuples (read1, read2)
  """
  ref_dict = dict([(r, n) for n, r in enumerate(fp.references)] + [('*', -1)])
  #                                                          unmapped with no contig
  ri = read_iter(fp, contig_q)
  singles = OrderedDict()
  while 1:
    if is_singles_mixer:
      try:
        read_str = singles_q.get_nowait()  # Any singles hanging about?
        if read_str is None: # One process says they are done with singles
          single_src_cnt -= 1
          if single_src_cnt == 0:  # Everyone says they've sent all their singles
            read = None
          else:
            continue  # At least one more source of singles about
        else:
          read = fromstring(read_str, ref_dict)
      except queue.Empty:
        read = next(ri, None)
        if read is None:
          time.sleep(0.01)  # We are out of contigs and we should keep checking the singles Q
          continue
    else:
      read = next(ri, None)

    if read is None:  # Out of reads from contigs and, if we are a mixer, out of reads from singles_q
      break

    key = read.qname[:20]  # Is this enough?
    if key not in singles:
      singles[key] = read
    else:  # You complete me
      yield (read, singles[key]) if read.is_read1 else (singles[key], read)
      del singles[key]

    if not is_singles_mixer:
      if len(singles) > max_singles:  # Flush earliest singles
        singles_q.put(singles.popitem(last=False).tostring(fp))

  # We need to send the remaining singles to the mixer
  if not is_singles_mixer:
    for read in singles.values():
      singles_q.put(read.tostring(fp))
    singles_q.put(None)  # No more singles from us


def worker(pipeline, bam_fname, result_q, contig_q,
           paired=False, singles_q=None, max_singles=1000,
           is_singles_mixer=False, single_src_cnt=None):
  """Given a pipeline, run it with reads from the given bam taken from contigs supplied
  over the contig_q.

  This expects the pipeline to yield one final result which it can then return.

  It expects the last element of pipeline to be a function that consumes a read iterator and returns a result.

  This is more flexible than you think, since the result can be an iterator, so this can be
  used to filter reads in parallel. See examples in the filter analysis tutorial


  :param pipeline:  A list of pipelines
  :param bam_fname: Source BAM file
  :param result_q:  The result is put here.
  :param contig_q:  messages are of the form (ref, True/False)
                    ref is the name of the contig
                    True/False indicates if eof should be set T/F
                    This controls whether we read to end of file including all the
                    unmapped reads. The caller figures out if this is that last
                    contig that sits just before that tail of unmapped reads at the end
                    of the BAM file
  :param paired:    Do we pair the reads before passing them to the pipeline?
  :param singles_q: messages are SAM strings of reads converted using tostring().
                    This is only used/relevant if paired=True because we use that to
                    collect the singles from all contigs and pair them up
                    Depending on whether this is the last

  :param max_singles:       When we have these many singles, start passing then to the
                            singles mixer

  :param is_singles_mixer:  Set True if this is also the "singles mixer" that
                            receives unpaired reads from other workers
  :param single_src_cnt:    How many sources of singles we have
                            This is
  :return:
  """
  if paired and singles_q is None:
    raise RuntimeError('Need singles_q to be defined if using paired reads')

  fp = pysam.AlignmentFile(bam_fname)
  if paired:
    t1 = paired_read_iter(fp, contig_q,
                          singles_q=singles_q, max_singles=max_singles,
                          is_singles_mixer=is_singles_mixer, single_src_cnt=single_src_cnt)
  else:
    t1 = unpaired_read_iter(fp, contig_q)

  sink = pipeline[-1]
  result_q.put(sink(cyt.pipe(t1, *pipeline[:-1])))


def scatter(pipeline, bam_fname, paired=False, ncpus=2, max_singles=1000):
  """Given a pipeline and a source bam file use multiprocessing to run the pipeline
  via multiple workers splitting up the work by contig

  python multiprocessing will be used for running the pipelines in parallel and care
  must be taken to ensure the individual pipeline nodes are parallelizable

  This expects the pipeline to yield one final result which it can then return.

  It expects the last element of pipeline to be a function that consumes a read iterator and returns a result.

  This is more flexible than you think, since the result can be an iterator, so this can be
  used to filter reads in parallel. See examples in the filter analysis tutorial

  :param bam_fname:
  :param pipeline:
  :param paired:  When run in parallel, paired vs unpaired pipelines work differently
                  So we have to tell scatter if we want to source paired or unpaired reads
  :param ncpus:
  :param max_singles:
  :return:
  """
  assert ncpus > 1, "ncpus = 1 can't use scatter!"

  result_q = Queue()
  contig_q = Queue()
  if paired:
    singles_q = Queue()
    is_mixer = [False] * (ncpus - 1) + [True]
  else:
    singles_q = False
    is_mixer = [False] * ncpus

  p_list = []
  for i in range(ncpus):
    p_list += [
      Process(target=worker,
              args=(pipeline, bam_fname, result_q, contig_q,
                    paired, singles_q, max_singles,
                    is_mixer[i], ncpus - 1))
    ]

  for p in p_list:
    p.start()

  _contigs = find_non_empty_contigs(bam_fname)
  contigs = [(c, False) for c in _contigs[:-1]] + [(_contigs[-1], True)]
  # This ensures that we read till EOF for the last contig and thereby fetch all of the trailing unmapped reads

  for contig in contigs:
    contig_q.put(contig)

  # Tell child processes to stop
  for i in range(ncpus):
    contig_q.put(None)

  for i in range(ncpus):
    yield result_q.get()

  # Orderly exit
  for p in p_list:
    p.join()


def find_non_empty_contigs(bam_fname):
  # Thanks to Güneş Bayir for suggesting a proper algorithm to pull the unmapped reads
  contigs = []
  fp = pysam.AlignmentFile(bam_fname)
  for ref in fp.references:
    for _ in fp.fetch(ref):
      contigs += [ref]
      break
  return contigs


def fromstring(s, ref_dict):
  """Inverse of pysam.AlignedSegment.tostring(): given a string, create an aligned segment

  :param s:
  :param ref_dict: ref_dict = dict([(r, n) for n, r in enumerate(fp.references)] + [('*', -1)])
  :return:
  """
  def _split(_s):
    qname, flag, rname, pos, \
    mapping_quality, cigarstring, \
    rnext, pnext, template_length, seq, qual, *_tg = _s.split('\t')

    flag = int(flag)
    rname = ref_dict[rname]  # dict must have '*': -1 entry too
    pos = int(pos)
    mapping_quality = int(mapping_quality)
    rnext = rname if rnext == '=' else ref_dict[rnext]
    pnext = int(pnext)
    template_length = int(template_length)

    return qname, flag, rname, pos, \
      mapping_quality, cigarstring, \
      rnext, pnext, template_length, seq, qual, _tg

  # So close, pysam.tostring, so close
  def _tags(_t):
    _tl = _t.split(':')
    if _tl[1] == 'i':
      _tl[2] = int(_tl[2])
    elif _tl[1] == 'f':
      _tl[2] = float(_tl[2])

    return _tl[0], _tl[2], _tl[1]

  r = pysam.AlignedSegment()
  r.qname, r.flag, r.rname, r.pos, \
  r.mapping_quality, r.cigarstring, \
  r.rnext, r.pnext, r.template_length, r.seq, r.qual, tags = _split(s)
  r.set_tags([_tags(t) for t in tags])
  return r