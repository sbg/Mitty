import logging

import click


import mitty.simulation.readgen as rgen


@click.group()
@click.version_option()
def cli():
  """A genomic data simulator for testing and debugging bio-informatics tools"""
  logging.basicConfig(level=logging.WARNING)
  pass


def print_qname(ctx, param, value):
  if not value or ctx.resilient_parsing:
    return
  click.echo(rgen.__qname_format_details__)
  ctx.exit()


@cli.command(short_help='Generate simulated reads.')
@click.argument('reffasta')
@click.argument('paramfile')
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
@click.option('--seed', default=7, help='Seed for RNG')
@click.option('--reference-reads-only', is_flag=True, help='If this is set, only reference reads will be generated')
@click.option('--sample-name', '-s', type=str, help='Name of sample to put in reads')
@click.option('--vcf', type=click.Path(exists=True), help='Path to vcf file representing sample')
@click.option('--bed', type=click.Path(exists=True), help='Bed file restricting regions reads are taken from')
def reads(reffasta, paramfile, sample_name, vcf, bed, reference_reads_only):
  """Generate simulated reads and write them to stdout"""
  pass


@cli.command(short_help='PE -> SE')
def unpair_fastq(fastq):
  """Rewrite qnames from interleaved PE FASTQ so file can be run as SE. Output to stdout"""
  pass


@cli.command(short_help='Interleaved PE 3> PE_1 4> PE_2')
@click.argument('fastq')
def uninterleave_fastq(fastq):
  """Given an interleaved paired-end FASTQ, split it out into two FASTQ files, written to streams 3 and 4"""
  pass


@cli.command()
def genomes():
  """Functions for sample/population simulations"""
  pass


