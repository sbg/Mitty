"""

file0 = A
file1 = B
file2 = C

ABC = present in all partitions
abc = absent from all partitions

"""
import os
import time
import logging

import pysam

from mitty.benchmarking.alignmentscore import ri, load_qname_sidecar, parse_qname, score_alignment_error, tag_alignment

logger = logging.getLogger(__name__)

MAX_ORIGINS = 8
labels = [['(A)', '(B)', '(C)', '(D)', '(E)', '(F)', '(G)', '(H)'],
          ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']]


def iterate_over_bams(bam_fp_l):
  """Give us an iterator that returns us a read (and which BAM it's from) until we run out of reads from every BAM

  :param bam_fp_l:
  :return:
  """
  bl = [bam_fp.fetch(until_eof=True) for bam_fp in bam_fp_l]
  while 1:
    rl = [next(bam_fp, None) for bam_fp in bl]
    if any(rl):
      for n, r in enumerate(rl):
        if r is None: continue
        if r.flag & 0x900 != 0: continue  # Not primary alignment
        yield n, r
    else:
      break


def main(bam_in_l, out_prefix, criterion, threshold, sidecar_fname=None):
  """

  :param bam_in_l:
  :param out_prefix:
  :param criterion: {'d_err', 'MQ', 'mapped', 'p_diff'}
  :param threshold:
  :param simulated:
  :param sidecar_fname
  :return:
  """
  assert len(bam_in_l) <= MAX_ORIGINS, "Can't do more than {} sets".format(MAX_ORIGINS)

  bam_fp_l = [pysam.AlignmentFile(bam_in) for bam_in in bam_in_l]
  long_qname_table = load_qname_sidecar(sidecar_fname) if sidecar_fname else None

  part_d = get_partition_description(out_prefix, len(bam_in_l))
  for p in part_d:
    p['filehandles'] = [pysam.AlignmentFile(p['filenames'][k] + '.unsorted.bam', 'wb', header=bam_fp_l[k].header)
                        for k in range(len(bam_fp_l))]
  scoring_fn = scoring_fn_dict.get(criterion)[0]

  incomplete_reads = {}
  cnt = -1
  t0 = time.time()
  for cnt, (n, r) in enumerate(iterate_over_bams(bam_fp_l)):
    if (cnt + 1) % 1000000 == 0:
      t1 = time.time()
      logger.debug('Processed {} reads ({} incomplete ({:0.2}%)) in {:.2f}s ({:.2f} t/s)'.format(
        cnt + 1, len(incomplete_reads), 100 * len(incomplete_reads) / (cnt + 1), t1 - t0, cnt / (t1 - t0)))

    ky = ('1' if r.is_read1 else '2') + r.qname

    if ky in incomplete_reads:
      ir = incomplete_reads[ky]
      ir[n] = r
      if all(ir):
        process_these_reads(part_d, incomplete_reads.pop(ky), scoring_fn, threshold, long_qname_table)
    else:
      ir = [None] * len(bam_fp_l)
      ir[n] = r
      incomplete_reads[ky] = ir

  t1 = time.time()
  logger.debug('Processed {} reads ({} incomplete ({:0.2}%)) in {:.2f}s ({:.2f} t/s)'.format(
    cnt + 1, len(incomplete_reads), 100 * len(incomplete_reads)/(cnt + 1), t1 - t0, cnt / (t1 - t0)))

  logger.debug('Closing output files')
  for p in part_d:
    for fp in p['filehandles']:
      fp.close()

  # Nice to get this written out before the time consuming sort and index stages
  with open('{}_summary.txt'.format(out_prefix), 'w') as fp:
    for p in part_d:
      fp.write('{}\t{}\n'.format(p['partition_label'], p['total']))

  logger.debug('Sorting and indexing output BAMs')
  for p in part_d:
    for fn in p['filenames']:
      logger.debug('Sort and index {}'.format(fn))
      pysam.sort('-m', '1G', '-o', fn + '.bam', fn + '.unsorted.bam')
      os.remove(fn + '.unsorted.bam')
      pysam.index(fn + '.bam')


def get_partition_description(file_prefix, n_bams):
  """Compute partitions and the files that correspond to each partition

  :param n_bams: number of original files 1, 2, ...
  :param simulated:
  :return:
  """
  def _part(_n):
    p_label = ''.join(get_partition_label(_n, n_bams)[::-1])
    return {
      'partition_label': p_label,
      'total': 0,
      'filenames': [
        '_'.join([file_prefix, p_label, labels[1][k]])
        for k in range(n_bams)
      ]
    }

  return [_part(n) for n in range(2 ** n_bams)]


def process_these_reads(part_d, rl, scoring_fn, threshold, long_qname_table):
  p = scoring_fn(rl, threshold, long_qname_table)
  if p is not None:
    # several of these scoring functions mutate reads in rl
    for bam_fp, r in zip(part_d[p]['filehandles'], rl):
      bam_fp.write(r)
    part_d[p]['total'] += 1


def get_partition_label(n, no):
  if no == 0:
    return []
  return [labels[n % 2][no - 1]] + get_partition_label(int(n / 2), no - 1)


# Partition functions ---------------------------------------------------------------
# The functions return a number between 0 and 2^n - 1 indicating the partition the read should
# go into. If the result is None, this read should be discarded (does not contribute to any counts)

def partition_by_d_err(rl, x, long_qname_table):
  """

  :param rl: a list of read objects [0, 1, 2, ...]
  :param x: threshold
  :return:
  """
  ri = parse_qname(rl[0].qname, long_qname_table)[0 if rl[0].is_read1 else 1]
  for r in rl:
    tag_alignment(r, ri)
  return independent_partitions(rl, lambda r: abs(r.get_tag('Xd')) < x)


def partition_by_d_err_strict(rl, x, long_qname_table):
  """

  :param rl: a list of read objects [0, 1, 2, ...]
  :param x: threshold
  :return:
  """
  ri = parse_qname(rl[0].qname, long_qname_table)[0 if rl[0].is_read1 else 1]
  for r in rl:
    tag_alignment(r, ri, strict=True)
  return independent_partitions(rl, lambda r: abs(r.get_tag('Xd')) < x)


def partition_by_MQ(rl, x, *args):
  """

  :param rl:
  :param x: MQ threshold
  :return:
  """
  return independent_partitions(rl, lambda r: r.mapping_quality < x)


def partition_by_mapped(rl, *args):
  """

  :param rl:
  :param args: rest are ignored
  :return:
  """
  return independent_partitions(rl, lambda r: not r.is_unmapped)


def partition_by_p_diff(rl, x, *args):
  """ Treating A as a gold standard, consider p_diff as the difference in alignment to A.
  Reads that are unmapped in A are discarded. Basically, modify the reads in place by
  filling out tag Xd then simply call partition_by_d_err for mapped reads. For unmapped
  reads return none

  :param rl: a list of read objects [0, 1, 2, ...]
  :param x: threshold
  :return:
  """
  r0 = rl[0]
  if r0.is_unmapped:
    logger.debug('{} is unmapped. Discarding for p_diff calcs'.format(r0.qname))
    return None

  # ['index', 'sample', 'chrom', 'cpy', 'strand', 'pos', 'cigar', 'special_cigar', 'v_list', 'md']
  r0_ri = ri(None, None, r0.reference_name, 0, 0, r0.pos, r0.cigarstring, None, None, None)
  for r in rl:
    tag_alignment(r, r0_ri)
  return independent_partitions(rl, lambda r: abs(r.get_tag('Xd')) < x)


pow2 = [2**n for n in range(MAX_ORIGINS)]


def independent_partitions(rl, scoring_function):
  """

  :param rl: a list of read objects [0, 1, 2, ...]
  :param scoring_function: a function that returns 0 or 1 to indicate if a read is out or in the
                           partition
  :return: a number between 0 and 2^n-1 where n = len(rl)
           indicating the partition
  """
  return sum([scoring_function(r) * pow2[n] for n, r in enumerate(rl[::-1])])


scoring_fn_dict = {
  'd_err': (partition_by_d_err,
            'Set membership is |d_err| < threshold\n'
            '(Only possible for simulated reads)'),
  'd_err_strict': (partition_by_d_err_strict,
            'Set membership is |d_err(strict)| < threshold\n'
            '(Only possible for simulated reads)'),
  'MQ': (partition_by_MQ,
         'Set membership is MQ < threshold'),
  'mapped': (partition_by_mapped,
             'Set membership is read.is_mapped'),
  'p_diff': (partition_by_p_diff,
             'Set membership is |p_diff| < threshold,\n'
             'where p_diff is computed by taking the reads of BAM A\n'
             'as the correct alignment, and ignoring unmapped reads in A')
}
