import logging
import os

import click


@click.group()
@click.version_option()
@click.option('-v', '--verbose', type=int, default=0)
def cli(verbose):
  """A genomic data simulator for testing and debugging bio-informatics tools"""
  logging.basicConfig(level=[
    logging.ERROR,
    logging.WARNING,
    logging.INFO,
    logging.DEBUG
  ][min(verbose, 3)])


@cli.command('filter-variants', short_help='Remove complex variants from VCF')
@click.argument('vcfin', type=click.Path(exists=True))
@click.argument('sample')
@click.argument('bed')
@click.argument('vcfout', type=click.Path())
def filter_vcf(vcfin, sample, bed, vcfout):
  """Subset VCF for given sample, apply BED file and filter out complex variants
   making it suitable to use for read generation"""
  import mitty.lib.vcfio as mvio
  mvio.prepare_variant_file(vcfin, sample, bed, vcfout)


@cli.command('gc-cov')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-b', '--block-len', type=int, default=10000, help='Block size for GC/cov computation')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def gc_cov(bam, fasta, pkl, block_len, threads):
  """Calculate GC content vs coverage from a BAM. Save in pickle file"""
  import mitty.empirical.gc as megc
  megc.process_bam_parallel(bam, fasta, pkl, block_len=block_len, threads=threads)


@cli.command('bq')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def sample_bq(bam, pkl, threads):
  """BQ distribution from BAM"""
  import mitty.empirical.bq as bbq
  bbq.process_bam_parallel(bam, pkl, threads=threads)


@cli.command('bam2illumina', short_help='Create read model from BAM file')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('pkl')
@click.argument('desc')
@click.option('--every', type=int, default=1, help='Sample every nth read')
@click.option('--min-mq', type=int, default=0, help='Discard reads with MQ less than this')
@click.option('-t', '--threads', type=int, default=2, help='Threads to use')
@click.option('--max-bp', type=int, default=300, help='Maximum length of read')
@click.option('--max-tlen', type=int, default=1000, help='Maximum size of insert')
def bam2illumina(bam, pkl, desc, every, min_mq, threads, max_bp, max_tlen):
  """Process BAM file and create an empirical model of template length and base quality
  distribution. The model is saved to a Python pickle file usable by Mitty read generation programs."""
  import mitty.empirical.bam2illumina as b2m
  b2m.process_bam_parallel(bam, pkl, model_description=desc, every=max(1, every), min_mq=min_mq,
                           threads=threads, max_bq=94, max_bp=max_bp, max_tlen=max_tlen)


@cli.command('list-read-models')
@click.option('-d', type=click.Path(exists=True), help='List models in this directory')
def list_read_models(d):
  """List read models"""
  import pkg_resources
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
      click.echo('\n----------\n{}:\n{}\n=========='.format(os.path.basename(mod_fname), mod_data['model_description']))
    except:
      logging.debug('Skipping {}. Not a read model file'.format(mod_fname))


@cli.command('describe-read-model')
@click.argument('modelfile')
@click.argument('figfile')
def describe_read_model(modelfile, figfile):
  """Plot panels describing this model"""
  read_module, model = get_read_model(modelfile)
  read_module.describe_model(os.path.basename(modelfile), model, figfile)


@cli.command()
def qname():
  """Display qname format"""
  from mitty.simulation.sequencing.writefastq import __qname_format_details__
  click.echo(__qname_format_details__)


def print_qname(ctx, param, value):

  if not value or ctx.resilient_parsing:
    return

  qname()
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
@click.argument('longqname', type=click.Path())
@click.option('--fastq2', type=click.Path())
@click.option('--threads', default=2)
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
def generate_reads(fasta, vcf, sample_name, bed, modelfile, coverage, seed, fastq1, longqname, fastq2, threads):
  """Generate simulated reads"""
  import mitty.simulation.readgenerate as reads

  read_module, model = get_read_model(modelfile)
  reads.process_multi_threaded(
    fasta, vcf, sample_name, bed, read_module, model, coverage,
    fastq1, longqname, fastq2, threads=threads, seed=seed)


@cli.command('corrupt-reads', short_help='Apply corruption model to FASTQ file of reads')
@click.argument('modelfile')
@click.argument('fastq1_in', type=click.Path(exists=True))
@click.argument('fastq1_out', type=click.Path())
@click.argument('sidecar_in', type=click.Path(exists=True))
@click.argument('sidecar_out', type=click.Path())
@click.argument('seed', type=int)
@click.option('--fastq2-in', type=click.Path(exists=True))
@click.option('--fastq2-out', type=click.Path())
@click.option('--threads', default=2)
def read_corruption(modelfile, fastq1_in, fastq1_out, sidecar_in, sidecar_out, seed, fastq2_in, fastq2_out, threads):
  """Apply corruption model to FASTQ file of reads"""
  import mitty.simulation.readcorrupt as rc

  read_module, read_model = get_read_model(modelfile)
  rc.multi_process(read_module, read_model, fastq1_in, fastq1_out, sidecar_in, sidecar_out,
                   fastq2_in, fastq2_out, processes=threads, seed=seed)


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
@cli.command('sampled-genome', short_help='Create a sample VCF from a population VCF')
@click.argument('vcfin', type=click.Path(exists=True))
@click.argument('vcfout', type=click.Path())
@click.option('--info-af', help='Population for AF')
@click.option('--af', type=float, help='fixed AF value if info-af not available')
@click.option('--bed', multiple=True, help='BED files to determine ploidy for each chromosome. See help')
def sampled_genome(vcfin, vcfout, info_af, af, bed):
  """The `sampled-genome` command can, given a VCF representing all variants in a population and their allele
frequencies, return a diploid sample VCF based on random sampling of the main list.

**Populations:**

  The program is passed an `info-af` parameter that looks for an INFO field of that name. The
  program interprets the value of that field as the alternative allele frequency and uses that to sample the
  respective variant. If you have a file, say like that from the 1000G project, that represents allele frequencies
  from multiple populations for each variant, you can simulate individuals from these different populations
  by selecting the appropriate tags. e.g. for the 1000G VCFs, `EUR_AF` for Europeans, `AMR_AF` for Americans etc.

  If this parameter is omitted the `af` parameter needs to be supplied. This sets a flat alternative allele
  frequency for all variants that is used for sampling.


**Ploidy and BED files:**

  The program accepts multiple BED files. A chromosome may appear zero or more times in each bed file.
  The ploidy of each simulated chromosome is the number of BED files in which it appears at least once.

  Say, for example, we have bed files that look like

\b
    1.bed:
    1 10  1000
    2 10  5000

\b
    2.bed:
    1 2000  3000

\b
    3.bed:
    2 1000 2000


If all three bed files are passed this would result in a simulated genome with two copies of both chrom 1 and chrom 2.
Due to the regions in the BED files, chrom1 would contain only HET variants while chrom2 may carry HOM variants in the
overlap region 1000-2000

If only 1.bed and 2.bed were passed the resulting simulated genome would have two copies of chrom 1 but one copy of
chrom 2"""
  pass


@cli.command('simulate-variants', short_help='Create a fully simulated VCF')
@click.argument('vcfout', type=click.File('w'))
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('sample')
@click.argument('bed', type=click.Path(exists=True))
@click.argument('seed', type=int)
@click.option('--p-het', default=0.6, type=float, help='Probability for heterozygous variants')
@click.option('--model', type=(str, float, int, int), multiple=True, help='<model type> <p> <min-size> <max-size>')
def simulate_variants(vcfout, fasta, sample, bed, seed, p_het, model):
  """Generates a VCF with simulated variants. The program carries three basic models for variant simulation
- SNPs, insertions and deletions and is invoked as follows:

\b
    mitty -v4 simulate-variants \
    - \  # Write the VCF to std out
    ~/Data/human_g1k_v37_decoy.fasta \
    mysample \ # The name of the sample to add to
    region.bed \
    7 \  # This is the random number generator seed
    --p-het 0.6 \   # The probability for heterozygous variants
    --model SNP 0.001 1 1 \   #  <model type> <p> <min-size> <max-size>
    --model INS 0.0001 10 100 \
    --model DEL 0.0001 10 100 | bgzip -c > sim.vcf.gz
  """
  import mitty.simulation.genome.simulatevariants as simvar
  simvar.main(fp_out=vcfout, fasta_fname=fasta, sample_name=sample, bed_fname=bed, seed=seed, p_het=p_het, models=model)


@cli.command('god-aligner', short_help='Create a perfect BAM from simulated FASTQs')
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('fastq1', type=click.Path(exists=True))
@click.argument('sidecar_in', type=click.Path(exists=True))
@click.argument('bam')
@click.option('--fastq2', type=click.Path(exists=True), help='If a paired-end FASTQ, second file goes here')
@click.option('--sample-name', default='S', help='If supplied, this is put into the BAM header')
@click.option('--platform-name', default='Illumina', help='If supplied, this is put into the BAM header')
@click.option('--max-templates', type=int, help='For debugging: quits after processing these many templates')
@click.option('--threads', default=2)
@click.option('--no-sort', is_flag=True, help='Leave the unsorted BAM fragments as is. Required if using an external tool to merge + sort + index')
def god_aligner(fasta, bam, sample_name, platform_name, fastq1, sidecar_in, fastq2, max_templates, threads, no_sort):
  """Given a FASTA.ann file and FASTQ made of simulated reads,
     construct a perfectly aligned BAM from them.

     A BAM produced by an aligner from the same FASTQ can be diff-d against the perfect BAM
     to check for alignment accuracy. (Also see the mitty filter-bam tool)

     The perfect BAM is also useful for testing variant callers by removing the aligner from the
     pipeline and reducing one moving part.

     Note: The program uses the fasta.ann file to construct the BAM header"""
  import mitty.benchmarking.god_aligner as god
  god.process_multi_threaded(
    fasta, bam, fastq1, sidecar_in, fastq2, threads, max_templates, platform_name, sample_name, not no_sort)


@cli.command('filter-bam', short_help='Refine alignment scoring')
def filter_bam():
  """Discard (simulated) reads from BAM whose alignment score passes a filter we set.

  This is meant to be run after running a BAM diff between an aligned BAM and a perfect BAM.
  The raw BAM diff will strictly flag any difference between the alignment and the perfect BAM
  but for our analyses it may be sufficient to be less strict.

  As an example, for some analysis we may not be concerned about alignments that have been soft-clipped
  but are otherwise correct. In other analyses we may want to include such soft-clipped reads because
  we suspect we are over-aggressive in soft-clipping.

  **Currently under development**
  """
  click.echo('Currently under development')


@cli.group('debug', short_help='Alignment and variant calling debugging tools')
def debug_tools():
  pass


@debug_tools.command('pr-by-size', short_help="Characterize TP, FN, FP and GT calls (and hence P/R) by variant size")
@click.argument('evcf', type=click.Path(exists=True))
@click.argument('out', type=click.Path())
@click.option('--region-label', help='Name of high confidence region if desired')
@click.option('--max-size', type=int, default=50, help='Maximum size of variant to consider')
@click.option('--title', help='Title for the plot')
@click.option('--fig-file', type=click.Path(), help='If supplied, plot will be saved here')
@click.option('--plot-bin-size', type=int, help='Bin size')
@click.option('--replot', is_flag=True,
              help='If supplied, instead of reprocessing the evcf, we expect "out" to exist, and load data from there')
def pr_by_size(evcf, out, region_label, max_size, title, fig_file, plot_bin_size, replot):
  import mitty.benchmarking.evcfbysize as ebs
  if not replot:
    data = ebs.main(evcf_fname=evcf, out_csv_fname=out,
                    max_size=max_size, high_confidence_region=region_label)
  else:
    data = ebs.np.loadtxt(out, skiprows=1, delimiter=',', dtype=[('TP', int), ('FN', int), ('GT', int), ('FP', int)])
  ebs.plot(data, fig_fname=fig_file, bin_size=plot_bin_size, title=title)


@debug_tools.command('variant-by-size', short_help="Characterize variant size distribution in a VCF")
@click.argument('vcf', type=click.Path(exists=True))
@click.argument('out', type=click.Path())
@click.option('--max-size', type=int, default=50, help='Maximum size of variant to consider')
@click.option('--title', help='Title for the plot')
@click.option('--fig-file', type=click.Path(), help='If supplied, plot will be saved here')
@click.option('--plot-bin-size', type=int, help='Bin size')
@click.option('--replot', is_flag=True,
              help='If supplied, instead of reprocessing the vcf, we expect "out" to exist, and load data from there')
def variant_by_size(vcf, out, max_size, title, fig_file, plot_bin_size, replot):
  import mitty.benchmarking.vsizedistrib as vsd
  if not replot:
    data = vsd.main(vcf_fname=vcf, max_size=max_size)
    vsd.np.savetxt(out, data, fmt='%d', delimiter=', ', header='SIZE')
  else:
    data = vsd.np.loadtxt(out, skiprows=1, delimiter=',', dtype=int)
  vsd.plot(data, fig_fname=fig_file, bin_size=plot_bin_size, title=title)


def partition_bam_choices():
  from mitty.benchmarking.partition_bams import scoring_fn_dict
  return scoring_fn_dict.keys()


def print_partion_bam_criteria(ctx, param, value):
  if not value or ctx.resilient_parsing:
    return
  from mitty.benchmarking.partition_bams import scoring_fn_dict
  for k, v in scoring_fn_dict.items():
    print('{}: {}\n'.format(k, v[1]))
  ctx.exit()


@debug_tools.command('partition-bams', short_help="Given two or more BAMs partition them into mutually exclusive sets")
@click.argument('outprefix')
@click.argument('criterion', type=click.Choice(partition_bam_choices()))
@click.option('--threshold', type=float, default=10)
@click.option('--sidecar_in', type=click.Path(exists=True))
@click.option('--bam', type=click.Path(exists=True), multiple=True, help='BAMs to partition')
@click.option('--criteria', is_flag=True, callback=print_partion_bam_criteria, expose_value=False, is_eager=True, help='Print documentation for criteria')
def partition_bams(outprefix, criterion, threshold, sidecar_in, bam):
  """
An example command line for this tool is:

\b
mitty -v4 debug partition-bams
  myderr \\
  d_err --threshold 10 \\
  --sidecar_in lq.txt --bam bwa_1.5.bam --bam bwa_10.bam --bam bwa_20.bam

This command line asks the tool to use |d_err| < 10 as the set membership function. We are passing it three BAM files
(the file names refer to the `-r` values we passed `bwa mem` (1.5, 10 and 20)) and `lq.txt` is the sidecar file carrying
the qnames > 254 characters (as described previously).

This tool produces a summary file `myderr_summary.txt` that looks like:

\b
(A)(B)(C) 22331
(A)(B)C   234
(A)B(C)   0
(A)BC     3
A(B)(C)   0
A(B)C     0
AB(C)     208
ABC       199126


In this nomenclature A is the set and (A) is the complement of this set. The set labels A, B, C ... (upto a maximum of 10)
refer to the BAM files in sequence, in this case 1.5, 10 and 20.

Thus, ABC means all the reads which have a |d_err| < 10 in all the three files. AB(C) means all the reads which have
a |d_err| < 10 in A and B but not C, and so on. A reader familiar with Venn diagrams is refered to the chart in the docs
for a translation of the three dimensional case to a three way Venn diagram. Higher dimensions are harder to visualize
as Venn diagrams.

The tool also produces a set of files following the naming convention:

\b
myderr_(A)(B)(C)_A.bam
myderr_(A)(B)(C)_B.bam
myderr_(A)(B)(C)_C.bam
myderr_(A)(B)C_A.bam
myderr_(A)(B)C_B.bam
myderr_(A)(B)C_C.bam
...

The first part of the name follows the convention outlined above. The trailing A, B, C refer to the orginal source BAM of
the reads. So `myderr_(A)(B)(C)_B.bam` carries reads from bam B that have |d_err| >= 10 in all the three BAMs.

The criteria the `partition-bam` tool can be run on can be obtained by passing it the `--criteria` option.
  """
  import mitty.benchmarking.partition_bams as pbm
  pbm.main(bam_in_l=bam, out_prefix=outprefix, criterion=criterion, threshold=threshold, sidecar_fname=sidecar_in)


@debug_tools.command('mq-plot', short_help='Create MQ plot from BAM')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('csv')
@click.argument('plot')
@click.option('--max-d', type=int, default=200, help='Range of d_err to consider')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--sample-name', help='If the FASTQ contains multiple samples, process reads from only this sample')
@click.option('--threads', default=2)
def mq_plot(bam, csv, plot, max_d, strict_scoring, sample_name, threads):
  """Expects the BAM to contain a tag called 'Xd' which carries an integer indicating
  how far off the alignment is."""
  import mitty.benchmarking.mq as mq
  mq_mat = mq.process_bam(bam, max_d, strict_scoring, sample_name, threads, csv)
  mq.plot_mq(mq_mat, plot)


@debug_tools.command('alignment-analysis-plot', short_help='Plot various alignment metrics from BAM')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('long_qname_file', type=click.Path(exists=True))
@click.argument('out', type=click.Path())
@click.option('--max-d', type=int, default=200, help='Range of d_err to consider')
@click.option('--max-size', type=int, default=50, help='Maximum size of variant to consider')
@click.option('--fig-file', type=click.Path(), help='If supplied, plot will be saved here')
@click.option('--plot-bin-size', type=int, help='Bin size')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--replot', is_flag=True,
              help='If supplied, instead of reprocessing the evcf, we expect "out" to exist, and load data from there')
@click.option('--processes', default=2, help='How many processes to use for computation')
def alignment_debug_plot(bam, long_qname_file, out, max_d, max_size, fig_file, plot_bin_size, strict_scoring, replot, processes):
  """Computes 3D matrix of alignment metrics (d_err, MQ, v_size) saves it to a numpy array file and produces a set
  of summary figures"""
  import mitty.benchmarking.xmv as xmv
  if not replot:
    xmv_mat = xmv.main(bam, sidecar_fname=long_qname_file,
                       max_xd=max_d, max_vlen=max_size, strict_scoring=strict_scoring,
                       processes=processes)
    xmv.save(xmv_mat, out)
  else:
    xmv_mat = xmv.np.load(out)

  xmv.plot_panel(xmv_mat, fig_fname=fig_file)


@cli.command('derr-plot', short_help='Create alignment error plot from scored BAM')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('longqname', type=click.Path(exists=True))
@click.argument('csv')
@click.argument('plot')
@click.option('--max-v', type=int, default=200, help='Range of variant sizes to consider (51 min)')
@click.option('--max-d', type=int, default=200, help='Range of d_err to consider')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--sample-name', help='If the FASTQ contains multiple samples, process reads from only this sample')
@click.option('--threads', default=2)
def derr_plot(bam, longqname, csv, plot, max_v, max_d, strict_scoring, sample_name, threads):
  """Expects the BAM to contain a tag called 'Xd' which carries an integer indicating
  how far off the alignment is."""
  import mitty.benchmarking.derr as derr
  derr_mat = derr.process_bam(bam, longqname, max_v, max_d, strict_scoring, sample_name, threads, csv)
  derr.plot_derr(derr_mat, plot)


@cli.command('filter-eval-vcf', short_help='Split out the FP and FN from an eval.vcf')
@click.argument('vcfin', type=click.Path(exists=True))
@click.argument('outprefix')
def filter_eval_vcf(vcfin, outprefix):
  """Subset VCF for given sample, apply BED file and filter out complex variants
   making it suitable to use for read generation"""
  import mitty.benchmarking.filterevalvcf as fev
  fev.extract_fp_fn(vcfin, outprefix)


def get_read_model(modelfile):
  """Return read module and model data given modelfile

  :param modelfile:
  :return:
  """
  import pickle
  import pkg_resources

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