import logging
import os
import json

import click


import mitty.lib.vcfio as mvio
import mitty.empirical.bq as bbq
import mitty.empirical.bq_fastq as bbqf
import mitty.empirical.gc as megc

import mitty.simulation.reads as reads
import mitty.simulation.illumina  # Hard coded for now, might use entry points like before to pip install models

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


@cli.command('bq-fastq')
@click.argument('fastq1', type=click.Path(exists=True))
@click.argument('fastq2', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def sample_bq_fastq(fastq1, fastq2, pkl, threads):
  """BQ distribution from BAM"""
  bbqf.base_quality(fastq1, fastq2, pkl, threads=threads)


@cli.command()
def qname():
  """Display qname format"""
  click.echo(reads.__qname_format_details__)


def print_qname(ctx, param, value):
  if not value or ctx.resilient_parsing:
    return
  click.echo(reads.__qname_format_details__)
  ctx.exit()


@cli.command('generate-reads', short_help='Generate simulated reads.')
@click.argument('modelfile')
@click.argument('fasta')
@click.argument('vcf')
@click.argument('sample_name')
@click.argument('bed')
@click.argument('seed', type=int)
@click.option('--fastq1')
@click.option('--fastq2')
@click.option('--threads', default=2)
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
def generate_reads(modelfile, fasta, vcf, sample_name, bed, seed, fastq1, fastq2, threads):
  """Generate simulated reads"""
  read_module = mitty.simulation.illumina
  # hard coded for now. In the future this will be read from the modelfile as before
  model_params = json.load(open(modelfile, 'r'))
  reads.process_multi_threaded(
    fasta, vcf, sample_name, bed, read_module, model_params,
    fastq1, fastq2, threads=threads, seed=seed)


@cli.command('corrupt-reads', short_help='Apply corruption model to FASTQ file of reads')
def read_corruption():
  pass


def read_bias():
  pass


@cli.command('unpair-fastq', short_help='PE -> SE')
def unpair_fastq(fastq):
  """Rewrite qnames from interleaved PE FASTQ so file can be run as SE. Output to stdout"""
  pass


@cli.command('uninterleave-fastq', short_help='Interleaved PE 3> PE_1 4> PE_2')
@click.argument('fastq')
def uninterleave_fastq(fastq):
  """Given an interleaved paired-end FASTQ, split it out into two FASTQ files, written to streams 3 and 4"""
  pass


@cli.command()
def genomes():
  """Functions for sample/population simulations"""
  pass


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
def mq_plot():
  """Given a BAM from simulated reads, construct an MQ plot"""
  pass