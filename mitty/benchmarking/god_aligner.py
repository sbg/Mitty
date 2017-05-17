import logging
import os
import sys
import time
from multiprocessing import Process, Queue

import pysam

from mitty.lib.cigars import cigarv2_v1
from mitty.simulation.readgenerate import DNA_complement
from mitty.simulation.sequencing.writefastq import load_qname_sidecar, parse_qname
from mitty.version import __version__

logger = logging.getLogger(__name__)

__process_stop_code__ = 'SETECASTRONOMY'


def construct_header(fasta_ann, rg_id, sample='S', platform='Illumina'):
  return {
    'HD': {'VN': '1.0'},
    'PG': [{'CL': ' '.join(sys.argv),
            'ID': 'mitty-god-aligner',
            'PN': 'god-aligner',
            'VN': __version__}],
    'RG': [{'ID': rg_id,
            'CN': 'MittyReadSimulator',
            'PL': platform,
            'SM': sample}],
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


def process_multi_threaded(
  fasta, bam_fname, fastq1, sidecar_fname, fastq2=None, threads=1, max_templates=None,
  platform='Illumina',
  sample_name='Seven',
  cigar_v2=True,
  do_not_index=False):
  """

  :param bam_fname:
  :param bam_hdr:
  :param fastq1:
  :param sidecar_fname: File containing just the long qnames
  :param fastq2:
  :param threads:
  :param max_templates:
  :param platform:
  :param sample_name:
  :param cigar_v2: If True, write out CIGARs in V2 format
  :param do_not_index: If True, the output BAMs will be collated into one bam, sorted and indexed
                       the N output BAMs created by the individual workers will be deleted at the end.
                       If False, the N output BAMs created by the individual workers will remain. This
                       option allows users to merge, sort and index the BAM fragments with their own tools
  :return:

  Note: The pysam sort invocation expects 1GB/thread to be available
  """
  long_qname_table = load_qname_sidecar(sidecar_fname)

  rg_id = 'rg{}'.format(hash(' '.join(sys.argv)))
  fasta_ann = fasta + '.ann'
  bam_hdr = construct_header(fasta_ann, rg_id=rg_id, sample=sample_name, platform=platform)

  # Start worker processes
  logger.debug('Starting {} processes'.format(threads))
  file_fragments = ['{}.{:04}.bam'.format(bam_fname, i) for i in range(threads)]

  in_queue = Queue()
  p_list = [Process(target=disciple,
                    args=(file_fragments[i], bam_hdr, rg_id, long_qname_table, cigar_v2,
                          in_queue))
            for i in range(threads)]
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

  merge_sorted_fragments(bam_fname, file_fragments, do_not_index=do_not_index)


def disciple(bam_fname, bam_hdr, rg_id, long_qname_table, cigar_v2, in_queue):
  """Create a BAM file from the FASTQ lines fed to it via in_queue

  :param bam_fname:
  :param bam_hdr:
  :param rg_id:
  :param long_qname_table:
  :param cigar_v2:
  :param in_queue:
  :return:
  """
  logger.debug('Writing to {} ...'.format(bam_fname))
  t0 = time.time()
  fp = pysam.AlignmentFile(bam_fname, 'wb', header=bam_hdr)
  ref_dict = {k['SN']: n for n, k in enumerate(bam_hdr['SQ'])}
  cnt = 0
  for cnt, (qname, read_data) in enumerate(iter(in_queue.get, __process_stop_code__)):
    write_perfect_reads(qname, rg_id, long_qname_table, ref_dict, read_data, cigar_v2, fp)
  fp.close()
  t1 = time.time()
  logger.debug('... {}: {} reads in {:0.2f}s ({:0.2f} t/s)'.format(bam_fname, cnt, t1 - t0, cnt/(t1 - t0)))

  logger.debug('Sorting {} -> {}'.format(bam_fname, bam_fname + '.sorted'))
  t0 = time.time()
  pysam.sort('-m', '1G', '-o', bam_fname + '.sorted', bam_fname)
  os.remove(bam_fname)
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('Shutting down thread for {}'.format(bam_fname))


def write_perfect_reads(qname, rg_id, long_qname_table, ref_dict, read_data, cigar_v2,
                        fp):
  """Given reads begining to a template, write out the perfect alignments to file

  :param qname:
  :param rg_id:
  :param long_qname_table:
  :param ref_dict: dict containing reference names mapped to ref_id
  :param read_data: [x1, x2, ... ] where xi is (seq, qual)
                     and, e.g. [x1, x2] constitute a pair, if the input is paired end
  :param cigar_v2: output CIGARs in V2 format
  :param fp: pysam file pointer to write out
  :return:
  """
  reads = [
    pysam.AlignedSegment()
    for _ in range(len(read_data))
  ]

  for n, (ri, rd, read) in enumerate(zip(parse_qname(qname, long_qname_table), read_data, reads)):
    read.qname = qname
    read.reference_id = ref_dict[ri.chrom]
    read.pos = ri.pos - 1
    read.cigarstring = ri.cigar if cigar_v2 else cigarv2_v1(ri.cigar)
    read.mapq = 60
    read.set_tag('RG', rg_id, value_type='Z')

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


def merge_sorted_fragments(bam_fname, file_fragments, do_not_index=False):
  logger.debug('Merging sorted BAM fragments ...')
  t0 = time.time()
  pysam.merge('-rpcf', bam_fname, *[f + '.sorted' for f in file_fragments])
  t1 = time.time()
  logger.debug('... {:0.2f}s'.format(t1 - t0))

  logger.debug('Removing fragments')
  for f in file_fragments:
   os.remove(f + '.sorted')

  if not do_not_index:
    logger.debug('BAM index ...')
    t0 = time.time()
    pysam.index(bam_fname, bam_fname + '.bai')
    t1 = time.time()
    logger.debug('... {:0.2f}s'.format(t1 - t0))