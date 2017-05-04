"""Contains the logic for handling read model invocation and read creation"""
import time
from multiprocessing import Process, Queue

import pysam
import numpy as np

import mitty.lib.vcfio as vio
from mitty.lib.sanitizeseq import sanitize
import mitty.simulation.rpc as rpc
from mitty.simulation.sequencing.writefastq import writer

import logging

logger = logging.getLogger(__name__)


SEED_MAX = (1 << 32) - 1  # Used for seeding rng
__process_stop_code__ = 'SETECASTRONOMY'
DNA_complement = str.maketrans('ATCGN', 'TAGCN')

# Will hold the mem-shared read only data
vcf_df = None


# The sequence of steps is to
# 1. Generate read templates based on BED file  -> save as global
# 2. Start up a bunch of workers and pass regions to them
# 3. Each worker creates the nodelist for the region and uses the relavant templates to create reads
# 4. Write out the reads as they are generated
#
# This code generates perfect reads. The reads can be piped (Tee-d) to an algorithm that implements
# corruption

# The following functions save data to global (Read only) variables so that we can use them with
# multi-processing

#TODO: Figure out behavior for SE model
#TODO: Add to readme the note that a single ended model will still produce two files, r2 will be empty
def process_multi_threaded(fasta_fname, vcf_fname, sample_name, bed_fname,
                           read_module, model, coverage,
                           fastq1_fname, sidecar_fname, fastq2_fname,
                           truncate_to=None,
                           unpair=False,
                           threads=2, seed=7):
  """

  :param fasta_fname:
  :param vcf_fname:
  :param sample_name:
  :param bed_fname:
  :param read_module:
  :param model:
  :param coverage:
  :param fastq1_fname:
  :param sidecar_fname:
  :param fastq2_fname:
  :param truncate_to: If set, shorten reads to this length
  :param unpair: If True, for models that produce paired-end reads, produce single end, instead
  :param threads:
  :param seed:
  :return:
  """
  if truncate_to is not None:
    model['mean_rlen'] = min(truncate_to, model['mean_rlen'])
    logger.debug('Truncating mean read length to: {}bp'.format(model['mean_rlen']))
  if unpair:
    logger.debug('Un-pairing reads')
    model['unpaired'] = True
    assert fastq2_fname is None, 'For unpaired reads, no FASTQ2 should be supplied'
  read_model = read_module.read_model_params(model, coverage)

  global vcf_df
  vcf_df = vio.load_variant_file(vcf_fname, sample_name, bed_fname)

  in_queue = Queue()
  out_queue = Queue(10000)  # This is the key to preventing workers from dying

  logger.debug('Starting {} workers'.format(threads))
  workers = []
  for worker_id in range(threads):
    p = Process(target=read_generating_worker,
                args=(worker_id, fasta_fname, sample_name, read_module, read_model, in_queue, out_queue))
    p.start()
    workers.append(p)

  logger.debug('Starting writer process')
  wr = Process(target=writer, args=(fastq1_fname, sidecar_fname, fastq2_fname, out_queue))
  wr.start()

  logger.debug('Loading up the work queue')
  for wd in get_data_for_workers(read_model, vcf_df, seed):
    in_queue.put(wd)

  logger.debug('Stopping workers')
  for _ in workers:
    in_queue.put(__process_stop_code__)
  for w in workers:
    w.join()
  logger.debug('All workers are done')

  logger.debug('Stopping writer')
  out_queue.put(__process_stop_code__)
  wr.join()


def get_data_for_workers(model, vcf, seed):
  """Given a vcf structure return an iterator that keeps returning us a region and copy to generate
  reads from until the desired coverage has been attained.

  :param model: Must contain a 'passes' field that determines how many times a region is passed to it
  :param vcf: ploidy aware variant structure as returned by vcfio. We only use the region and ploidy information
              to create and return a generator for regions. We generate as many regions as needed to achieve the
              given coverage, taking ploidy of the region into account
  :param seed: This is used to seed an RNG that generates rng seeds passed to all the workers
  :return: a generator that returns worker data as a dict
           {
             'region_idx': r_idx,  # Index into VCF structure
             'region_cpy': cpy,  # which copy of the region to use
             'rng_seed': s   # seed for the RNGs in the worker
           }

           access the variants as vcf_df[r_idx]['v'][cpy]
  """
  seed_rng = np.random.RandomState(seed)
  shuffle_seed = seed_rng.randint(SEED_MAX)
  region_list = [
    {'region_idx': idx, 'region_cpy': cpy, 'rng_seed': seed_rng.randint(SEED_MAX)}
    for idx, v in enumerate(vcf)
    for cpy in range(len(v['v']))
    for _ in range(model['passes'])
  ]
  logger.debug('{} passes will be made'.format(len(region_list)))
  np.random.RandomState(shuffle_seed).shuffle(region_list)

  for region in region_list:
    yield region


# def read_generating_worker(worker_id, fasta_fname, sample_name, read_module, read_model, in_queue, out_queue):
#   for ps, wd in enumerate(iter(in_queue.get, __process_stop_code__)):
#     out_queue.put([str(ps), str(wd)])


def read_generating_worker(worker_id, fasta_fname, sample_name, read_module, read_model, in_queue, out_queue):
  """This worker will be given a fasta_fname, and region information. It is to generate the node_list,
  the read locations and then, finally, the reads themselves. The reads - in FASTQ format - are returned
  to the parent thread for writing.

  :param worker_id: Just a serial that helps us number reads uniquely
  :param fasta_fname:
  :param sample_name:  (This is just passed to the read qname)
  :param model:
  :param region_idx:
  :param in_queue: will be sending work information as (region_idx, cpy)
  :param out_queue: We'll be sending strings representing each template in FASTQ format
  :return:
  """
  fasta = pysam.FastaFile(fasta_fname)
  total_cnt, t00 = 0, time.time()
  for ps, wd in enumerate(iter(in_queue.get, __process_stop_code__)):
    r_idx, cpy, rng_seed = wd['region_idx'], wd['region_cpy'], wd['rng_seed']
    region = vcf_df[r_idx]['region']
    ref_seq = sanitize(fasta.fetch(reference=region[0], start=region[1], end=region[2]))

    # The structure needed to generate read sequences, pos and CIGAR strings
    # + 1 because BED is 0-indexed and our convention is 1-indexed as is what is displayed in genome browsers
    node_list = rpc.create_node_list(ref_seq, ref_start_pos=region[1] + 1, vl=vcf_df[r_idx]['v'][cpy])

    p_min, p_max = node_list[0].ps, node_list[-1].ps + node_list[-1].oplen  # We never end with a deletion
    r_info_l = read_module.generate_reads(read_model, p_min, p_max, rng_seed)

    qname_serial_stub = '{}:{}:{}'.format(sample_name, worker_id, ps)
    t0 = time.time()
    this_cnt = 0
    for template in zip(
      *[zip(*([r_info[k] for k in ['file_order', 'strand', 'pos', 'len']] +
                rpc.get_begin_end_nodes(r_info['pos'], r_info['len'], node_list))) for r_info in r_info_l]):
      reads = [None] * len(template)
      for fo, s, p, l, ns, ne in template:
        pos, cigar, v_list, seq = rpc.generate_read(p, l, ns, ne, node_list)
        if seq.count('N') > 2: break
        # This break combined with the else clause skips those read pairs where
        # at least one read has too many 'N's
        if s == 1:
          seq = seq.translate(DNA_complement)[::-1]
        reads[fo] = (s, pos, cigar, v_list, '', seq, 'I' * len(seq))
      else:
        this_cnt += 1
        # (
        #   None,  - We want the writer to put a unique stamp on each template
        #   sample_name,
        #   chrom,
        #   copy,
        #   (
        #     (strand, pos, cigar, (v1,v2,...), MD, seq, qual)
        #     ...  [repeated as for as many reads in this template]
        #   )
        # )
        out_queue.put((None, sample_name, region[0], cpy, reads))

    t1 = time.time()
    logger.debug('Worker {} ({}): {} templates in {:0.2f}s ({:0.2f} t/s)'.format(worker_id, region, this_cnt, t1 - t0, this_cnt/(t1 - t0)))
    total_cnt += this_cnt

  t11 = time.time()
  logger.debug('Worker {} finished: {} templates in {:0.2f}s ({:0.2f} t/s)'.format(
    worker_id, total_cnt, t11 - t00, total_cnt / (t11 - t00)))


def validate_templates_from_read_model(tplt):
  """A series of tests to make sure that the returned template locations are in a format suitable for our use

  :param tplt:
  :return:
  """
  if not isinstance(tplt, list):
    raise RuntimeError('Read model should return a list')

  if not isinstance(tplt[0], dict):
    raise RuntimeError('Read model should return a list of dicts')

  for k in ['reads', 'region']:
    if k not in tplt[0]:
      raise RuntimeError('Read model dict should contain {}'.format(k))

  if not isinstance(tplt[0]['reads'], np.recarray):
    raise RuntimeError('Reads should be in the form of a numpy recarray')