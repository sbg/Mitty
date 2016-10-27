import logging
import pkg_resources
import os
import json

import click


import mitty.lib.vcfio as mvio
import mitty.empirical.bq as bbq
import mitty.empirical.gc as megc

import mitty.benchmarking.god_aligner as god


@click.group()
@click.version_option()
@click.option('-v', '--verbose', type=int, default=0)
def cli(verbose):
  """A genomic data simulator for testing and debugging bio-informatics tools"""
  if verbose == 1:
    logging.basicConfig(level=logging.ERROR)
  elif verbose == 2:
    logging.basicConfig(level=logging.WARNING)
  elif verbose == 3:
    logging.basicConfig(level=logging.INFO)
  elif verbose >= 4:
    logging.basicConfig(level=logging.DEBUG)


@cli.command('filter-variants', short_help='Remove complex variants from VCF')
@click.argument('vcfin', type=click.Path(exists=True))
@click.argument('sample')
@click.argument('bed')
@click.argument('vcfout', type=click.Path())
def filter_vcf(vcfin, sample, bed, vcfout):
  """Subset VCF for given sample, apply BED file and filter out complex variants
   making it suitable to use for read generation"""
  mvio.prepare_variant_file(vcfin, sample, bed, vcfout)


@cli.command('gc-cov')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-b', '--block-len', type=int, default=10000, help='Block size for GC/cov computation')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def gc_cov(bam, fasta, pkl, block_len, threads):
  """Calculate GC content vs coverage from a BAM. Save in pickle file"""
  megc.process_bam_parallel(bam, fasta, pkl, block_len=block_len, threads=threads)


@cli.command('bq')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def sample_bq(bam, pkl, threads):
  """BQ distribution from BAM"""
  bbq.process_bam_parallel(bam, pkl, threads=threads)


@cli.command('bam2illumina', short_help='Create read model from BAM file')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('pkl')
@click.argument('desc')
@click.option('--every', type=int, default=1, help='Sample every nth read')
@click.option('-t', '--threads', type=int, default=2, help='Threads to use')
@click.option('--max-bp', type=int, default=300, help='Maximum length of read')
@click.option('--max-tlen', type=int, default=1000, help='Maximum size of insert')
def bam2illumina(bam, pkl, desc, every, threads, max_bp, max_tlen):
  """Create read model from BAM file"""
  import mitty.empirical.bam2illumina as b2m
  b2m.process_bam_parallel(bam, pkl, model_description=desc, every=every, threads=threads, max_bq=94, max_bp=max_bp, max_tlen=max_tlen)


@cli.command('list-read-models')
@click.option('-d', type=click.Path(exists=True), help='List models in this directory')
def list_read_models(d):
  """List read models"""
  import pickle
  import glob

  if d is None:
    if not pkg_resources.resource_isdir(__name__, 'data/readmodels/'):
      logging.error('Read model directory not found in distribution!')
      raise FileNotFoundError
    model_list = [
      pkg_resources.resource_filename(__name__, 'data/readmodels/' + model)
      for model in pkg_resources.resource_listdir(__name__, 'data/readmodels/')]
  else:
    model_list = glob.glob(os.path.join(d, '*'))

  for mod_fname in model_list:
    try:
      mod_data = pickle.load(open(mod_fname, 'rb'))
      click.echo('{}: {}'.format(os.path.basename(mod_fname), mod_data['model_description']))
    except:
      logging.debug('Skipping {}. Not a read model file'.format(mod_fname))


@cli.command('describe-read-model')
@click.argument('modelfile')
@click.argument('figfile')
def describe_read_model(modelfile, figfile):
  """Plot panels describing this model"""
  read_module, model = get_read_model(modelfile)
  read_module.describe_model(model, figfile)


@cli.command()
def qname():
  """Display qname format"""
  import mitty.simulation.readgenerate as reads
  click.echo(reads.__qname_format_details__)


def print_qname(ctx, param, value):
  import mitty.simulation.readgenerate as reads

  if not value or ctx.resilient_parsing:
    return
  click.echo(reads.__qname_format_details__)
  ctx.exit()


@cli.command('generate-reads', short_help='Generate simulated reads.')
@click.argument('fasta')
@click.argument('vcf')
@click.argument('sample_name')
@click.argument('bed')
@click.argument('modelfile')
@click.argument('coverage', type=float)
@click.argument('seed', type=int)
@click.argument('fastq1', type=click.Path())
@click.option('--fastq2', type=click.Path())
@click.option('--threads', default=2)
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
def generate_reads(fasta, vcf, sample_name, bed, modelfile, coverage, seed, fastq1, fastq2, threads):
  """Generate simulated reads"""
  import mitty.simulation.readgenerate as reads

  read_module, model = get_read_model(modelfile)
  reads.process_multi_threaded(
    fasta, vcf, sample_name, bed, read_module, model, coverage,
    fastq1, fastq2, threads=threads, seed=seed)


@cli.command('corrupt-reads', short_help='Apply corruption model to FASTQ file of reads')
@click.argument('modelfile')
@click.argument('fastq1_in', type=click.Path(exists=True))
@click.argument('fastq1_out', type=click.Path())
@click.argument('seed', type=int)
@click.option('--fastq2-in', type=click.Path(exists=True))
@click.option('--fastq2-out', type=click.Path())
@click.option('--threads', default=2)
def read_corruption(modelfile, fastq1_in, fastq1_out, seed, fastq2_in, fastq2_out, threads):
  """Apply corruption model to FASTQ file of reads"""
  import mitty.simulation.readcorrupt as rc

  read_module, read_model = get_read_model(modelfile)
  rc.multi_process(read_module, read_model, fastq1_in, fastq1_out, fastq2_in, fastq2_out, processes=threads, seed=seed)


# def read_bias():
#   pass
#
#
# @cli.command('unpair-fastq', short_help='PE -> SE')
# def unpair_fastq(fastq):
#   """Rewrite qnames from interleaved PE FASTQ so file can be run as SE. Output to stdout"""
#   pass
#
#
# @cli.command('uninterleave-fastq', short_help='Interleaved PE 3> PE_1 4> PE_2')
# @click.argument('fastq')
# def uninterleave_fastq(fastq):
#   """Given an interleaved paired-end FASTQ, split it out into two FASTQ files, written to streams 3 and 4"""
#   pass
#
#
# @cli.command()
# def genomes():
#   """Functions for sample/population simulations"""
#   pass


@cli.command('god-aligner', short_help='Create a perfect BAM from simulated FASTQs')
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('fastq1', type=click.Path(exists=True))
@click.argument('bam')
@click.option('--fastq2', type=click.Path(exists=True))
@click.option('--sample-name')
@click.option('--max-templates', type=int)
@click.option('--threads', default=2)
def god_aligner(fasta, bam, sample_name, fastq1, fastq2, max_templates, threads):
  """Given a FASTA.ann file and FASTQ made of simulated reads,
     construct a perfectly aligned "god" BAM from them.

     This is useful for testing variant callers.

     The program uses the fasta.ann file to construct the BAM header"""
  god.process_multi_threaded(
    fasta, bam, fastq1, fastq2, threads, max_templates, sample_name)


@cli.command('mq-plot', short_help='Create MQ plot from BAM')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('csv')
@click.argument('plot')
@click.option('--max-d', type=int, default=200, help='Range of d_err to consider')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--threads', default=2)
def mq_plot(bam, csv, plot, max_d, strict_scoring, threads):
  """Given a BAM from simulated reads, construct an MQ plot"""
  import mitty.benchmarking.mq as mq
  mq_mat = mq.process_bam(bam, max_d, strict_scoring, threads, csv)
  mq.plot_mq(mq_mat, plot)


@cli.command('derr-plot', short_help='Create alignment error plot from BAM')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('csv')
@click.argument('plot')
@click.option('--max-v', type=int, default=200, help='Range of variant sizes to consider (51 min)')
@click.option('--max-d', type=int, default=200, help='Range of d_err to consider')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--threads', default=2)
def mq_plot(bam, csv, plot, max_v, max_d, strict_scoring, threads):
  """Given a BAM from simulated reads, show pattern of alignment errors"""
  import mitty.benchmarking.derr as derr
  derr_mat = derr.process_bam(bam, max_v, max_d, strict_scoring, threads, csv)
  derr.plot_derr(derr_mat, plot)


def get_read_model(modelfile):
  """Return read module and model data given modelfile

  :param modelfile:
  :return:
  """
  import pickle

  import mitty.simulation.illumina  # Hard coded for now, might use entry points like before to pip install models

  if pkg_resources.resource_exists(__name__, 'data/readmodels/' + modelfile):
    logging.debug('Found model {} in builtins'.format(modelfile))
    mod_fname = pkg_resources.resource_filename(__name__, 'data/readmodels/' + modelfile)
  else:
    logging.debug('Treating {} as literal path to model file'.format(modelfile))
    mod_fname = modelfile  # treat this as a literal file path

  model = pickle.load(open(mod_fname, 'rb'))
  read_module = {
    'illumina': mitty.simulation.illumina
  }.get(model['model_class'])

  return read_module, model