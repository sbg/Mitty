import logging
import os

import click


import mitty.lib.vcfio as mvio
import mitty.simulation.reads as reads
import mitty.empirical.bq as bbq
import mitty.empirical.bq_fastq as bbqf
import mitty.empirical.gc as megc


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
@click.argument('reffasta')
@click.argument('paramfile')
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
@click.option('--seed', default=7, help='Seed for RNG')
@click.option('--reference-reads-only', is_flag=True, help='If this is set, only reference reads will be generated')
@click.option('--sample-name', '-s', type=str, help='Name of sample to put in reads')
@click.option('--vcf', type=click.Path(exists=True), help='Path to vcf file representing sample')
@click.option('--bed', type=click.Path(exists=True), help='Bed file restricting regions reads are taken from')
def generate_reads(reffasta, paramfile, sample_name, vcf, bed, reference_reads_only):
  """Generate simulated reads and write them to stdout"""
  pass


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