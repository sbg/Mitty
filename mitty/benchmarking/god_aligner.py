import sys
import time
import base64
import logging
from multiprocessing import Process, Queue
import os

import pysam

from mitty.simulation.reads import parse_qname, DNA_complement
from mitty.version import __version__


logger = logging.getLogger(__name__)

__process_stop_code__ = 'SETECASTRONOMY'


def construct_header(fasta_ann, rg_id, sample='S'):
  return {
    'HD': {'VN': '1.0'},
    'PG': [{'CL': ' '.join(sys.argv),
            'ID': 'mitty-god-aligner',
            'PN': 'god-aligner',
            'VN': __version__}],
    'RG': [{'ID': rg_id, 'SM': sample}],
    'SQ': parse_ann(fasta_ann)
  }


def parse_ann(fn):
  """Given a fasta.ann file name parse it."""
  logger.debug('Parsing {} for sequence header information'.format(fn))
  ln = open(fn, 'r').readlines()[1:]  # Don't need the first line
  return [
    {
      'SN': ln[n].split()[1],
      'LN': int(ln[n + 1].split()[1])
    }
  for n in range(0, len(ln), 2)]


def process_multi_threaded(fasta, bam_fname, fastq1, fastq2=None, threads=1, max_templates=None, sample_name='Seven'):
  """

  :param bam_fname:
  :param bam_hdr:
  :param fastq1:
  :param fastq2:
  :param threads:
  :param max_templates:
  :return:
  """
  rg_id = base64.b64encode(' '.join(sys.argv).encode('ascii'))
  fasta_ann = fasta + '.ann'
  bam_hdr = construct_header(fasta_ann, rg_id=rg_id, sample=sample_name)

  # Start worker processes
  logger.debug('Starting {} processes'.format(threads))
  file_fragments = ['{}.{:04}'.format(bam_fname, i) for i in range(threads)]

  in_queue = Queue()
  p_list = [Process(target=disciple, args=(file_fragments[i], bam_hdr, in_queue)) for i in range(threads)]
  for p in p_list:
    p.start()

  t0 = time.time()

  # Burn through file
  logger.debug('Starting to read FASTQ file')
  fastq_l = [pysam.FastxFile(fastq1)]
  if fastq2 is not None: fastq_l += [pysam.FastxFile(fastq2)]

  cnt = 0
  for cnt, reads in enumerate(zip(*fastq_l)):
    # qname, [(seq, qual) ... ]
    in_queue.put((
        reads[0].name,
        [(r.sequence, r.quality) for r in reads]))
    if max_templates is not None and cnt >= max_templates:
      break
    if cnt % 100000 == 0:
      logger.debug('Read {} templates'.format(cnt))

  # Tell child processes to stop
  logger.debug('Stopping child processes')
  for i in range(threads):
    in_queue.put(__process_stop_code__)

  # Wait for them to finish
  for p in p_list:
    p.join()

  t1 = time.time()
  logger.debug('Processed {} templates in {:0.2f}s ({:0.2f} t/s)'.format(cnt, t1 - t0, cnt/(t1 - t0)))

  logger.debug('Combining BAM fragments ...')
  t0 = time.time()
  pysam.cat('-o', bam_fname + '.unsorted', *file_fragments)
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('BAM sort ...')
  t0 = time.time()
  pysam.sort('-o', bam_fname, bam_fname + '.unsorted')
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('BAM index ...')
  t0 = time.time()
  pysam.index(bam_fname)
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('Removing intermediate files')
  for f in file_fragments:
    os.remove(f)
  os.remove(bam_fname + '.unsorted')


def disciple(bam_fname, bam_hdr, in_queue):
  """Create a BAM file from the FASTQ lines fed to it via in_queue

  :param bam_fname:
  :param bam_hdr:
  :param in_queue:
  :return:
  """
  logger.debug('Writing to {}'.format(bam_fname))
  fp = pysam.AlignmentFile(bam_fname, 'wb', header=bam_hdr)
  ref_dict = {k['SN']: n for n, k in enumerate(bam_hdr['SQ'])}
  for qname, read_data in iter(in_queue.get, __process_stop_code__):
    write_perfect_reads(qname, ref_dict, read_data, fp)
  fp.close()
  logger.debug('Shutting down thread for {}'.format(bam_fname))


def write_perfect_reads(qname, ref_dict, read_data, fp):
  """Given reads begining to a template, write out the perfect alignments to file

  :param qname:
  :param read_data: [x1, x2, ... ] where xi is (seq, qual)
                     and, e.g. [x1, x2] constitute a pair, if the input is paired end
  :param fp: pysam file pointer to write out
  :return:
  """
  reads = [
    pysam.AlignedSegment()
    for _ in range(len(read_data))
  ]

  for n, (ri, rd, read) in enumerate(zip(parse_qname(qname), read_data, reads)):
    read.qname = qname
    read.reference_id = ref_dict[ri.chrom]
    read.pos = ri.pos - 1
    read.cigarstring = ri.cigar
    read.mapq = 60
    # TODO: ask aligner people what the graph cigar tag is
    # TODO: Set this as unmapped?
    if ri.strand:
      read.is_reverse = 1
      read.seq = rd[0].translate(DNA_complement)[::-1]
      read.qual = rd[1][::-1]
    else:
      read.is_reverse = 0
      read.seq = rd[0]
      read.qual = rd[1]

  # Some special stuff we have to do for paired reads
  if len(reads) == 2:
    for n, r in enumerate(reads):
      r.is_paired = True
      r.is_proper_pair = True
      r.is_read1 = (n == 0)
      r.is_read2 = (n == 1)
      r.pnext = reads[1 - n].pos
      r.rnext = reads[1 - n].reference_id

  for r in reads:
    fp.write(r)