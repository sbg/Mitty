"""Given a FASTQ (or BAM) and a long-qnames file trancate the qnames to
whatever length we want and push the too-long qnames to the overflow
file. Also convert any old style qnames (not ending in a *) to new style
ones with a proper termination character"""
import time
import logging

import pysam

from mitty.benchmarking.alignmentscore import load_qname_sidecar

logger = logging.getLogger(__name__)


def main(mainfile_in, sidecar_in, mainfile_out, sidecar_out, truncate_to=240, file_type=None):
  """

  :param mainfile_in:
  :param sidecar_in:
  :param mainfile_out:
  :param sidecar_out:
  :param truncate_to:
  :param file_type: If supplied ("BAM" or "FASTQ") then we will not auto detect
  :return:
  """
  ft = {
    'BAM': 1,
    'FASTQ': 0
  }.get(file_type or auto_detect(mainfile_in))
  fp_in = pysam.AlignmentFile(mainfile_in, mode='rb') if ft else pysam.FastxFile(mainfile_in)
  fp_out = pysam.AlignmentFile(mainfile_out, mode='wb', header=fp_in.header) if ft else open(mainfile_out, 'w')
  side_car_fp = open(sidecar_out, 'w')

  logger.debug('Starting conversion ...')
  long_qname_table = load_qname_sidecar(sidecar_in)
  cnt, t0 = 0, time.time()
  for cnt, r in enumerate(fp_in):
    qname = r.qname if ft else r.name  # Thanks pysam for the inconsistent naming. What's a few CPU cycles between friends?
    qname = long_qname_table.get(qname.split('|', 1)[0], qname)  # Don't pass the wrong side car file, you won't know what hit you
    qname = qname[:-1] + '*' # Older qnames had "|" instead of "*". "*" is more unambiguous as a termination character
    if len(qname) > truncate_to:
      side_car_fp.write('@' + qname + '\n')
      qname = qname[:truncate_to]
    if ft:
      r.qname = qname
      fp_out.write(r)
    else:
      fp_out.write('@{}\n{}\n+\n{}\n'.format(qname, r.sequence, r.quality))

    if cnt % 100000 == 99999:
      t1 = time.time()
      logger.debug('Processed {} reads in {:0.2f}s ({:0.2f} r/s)'.format(cnt + 1, t1 - t0, (cnt + 1)/(t1 - t0)))

  t1 = time.time()
  logger.debug('Processed {} reads in {:0.2f}s ({:0.2f} r/s)'.format(cnt + 1, t1 - t0, (cnt + 1) / (t1 - t0)))


def auto_detect(fname):
  return 'BAM' if fname.upper().endswith('BAM') else 'FASTQ'