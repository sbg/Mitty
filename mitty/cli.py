import logging
import os

import click


logger = logging.getLogger(__name__)


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

  import pkg_resources  # part of setuptools
  version = pkg_resources.require("Mitty")[0].version
  logger.debug('Mitty version {}'.format(version))


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


@cli.group('create-read-model', short_help='Create read model from different sources')
def create_read_model():
  pass


@create_read_model.command('bam2illumina', short_help='Create read model from BAM file')
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


@create_read_model.command('synth-illumina', short_help='Create fully synthetic Illumina-like read model')
@click.argument('pkl')
@click.option('--read-length', type=int, default=100)
@click.option('--mean-template-length', type=int, default=500)
@click.option('--std-template-length', type=int, default=50)
@click.option('--bq0', type=int, default=30, help='Maximum BQ')
@click.option('--k', type=float, default=200, help='Steepness of BQ curve. Larger => more steep')
@click.option('--sigma', type=float, default=10, help='Spread of BQ about mean value')
@click.option('--comment', help='Additional comments')
@click.option('--max-tlen', type=int, default=1000, help='Maximum template length we will generate')
def synth_read_model(pkl, read_length, mean_template_length, std_template_length, max_tlen, bq0, k, sigma, comment):
  """This generates a synthetic read model based on the parameters we pass it"""
  import mitty.simulation.sequencing.syntheticsequencer as synth
  synth.create_model(
    pkl,
    read_length=read_length,
    mean_template_length=mean_template_length, std_template_length=std_template_length, max_tlen=max_tlen,
    bq0=bq0, k=k, sigma=sigma, comment=comment)


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
@click.option('--truncate-to', type=int, help='Truncate all reads to these many bp (If set)')
@click.option('--unpair', is_flag=True, help='unpair reads for models that normally produce paired reads')
@click.option('--threads', default=2)
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
def generate_reads(fasta, vcf, sample_name, bed, modelfile,
                   coverage, seed,
                   fastq1, longqname, fastq2,
                   truncate_to, unpair,
                   threads):
  """Generate simulated reads"""
  import mitty.simulation.readgenerate as reads

  read_module, model = get_read_model(modelfile)
  reads.process_multi_threaded(
    fasta, vcf, sample_name, bed, read_module, model, coverage,
    fastq1, longqname, fastq2,
    truncate_to=truncate_to,
    unpair=unpair,
    threads=threads, seed=seed)


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


def print_variant_model_list(ctx, param, value):
  import mitty.simulation.genome.simulatevariants as simvar
  if not value or ctx.resilient_parsing:
    return

  for k in sorted(simvar.model_dispatch.keys()):
    print('{}:\n------'.format(k))
    print(simvar.model_dispatch[k].__doc__)
    print()
  ctx.exit()


@cli.command('simulate-variants', short_help='Create a fully simulated VCF')
@click.argument('vcfout', type=click.File('w'))
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('sample')
@click.argument('bed', type=click.Path(exists=True))
@click.argument('seed', type=int)
@click.option('--p-het', default=0.6, type=float, help='Probability for heterozygous variants')
@click.option('--model', type=(str, float, int, int), multiple=True, help='<model type> <p> <min-size> <max-size>')
@click.option('--list-models', is_flag=True, callback=print_variant_model_list, expose_value=False, is_eager=True, help='Print list of variant models')
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
@click.option('--cigar-v2', is_flag=True, help='Write out CIGARs in V2 for')
@click.option('--max-templates', type=int, help='For debugging: quits after processing these many templates')
@click.option('--threads', default=2)
@click.option('--do-not-index', is_flag=True, help='Leave the unsorted BAM fragments as is. Required if using an external tool to merge + sort + index')
def god_aligner(fasta, bam, sample_name, platform_name, fastq1, sidecar_in, fastq2,
                cigar_v2,
                max_templates,
                threads,
                do_not_index):
  """Given a FASTA.ann file and FASTQ made of simulated reads,
     construct a perfectly aligned BAM from them.

     A BAM produced by an aligner from the same FASTQ can be diff-d against the perfect BAM
     to check for alignment accuracy. (Also see the mitty filter-bam tool)

     The perfect BAM is also useful for testing variant callers by removing the aligner from the
     pipeline and reducing one moving part.

     Note: The program uses the fasta.ann file to construct the BAM header"""
  import mitty.benchmarking.god_aligner as god
  god.process_multi_threaded(
    fasta, bam, fastq1, sidecar_in, fastq2, threads, max_templates, platform_name, sample_name,
    cigar_v2=cigar_v2,
    do_not_index=do_not_index)


@cli.group('debug', short_help='Alignment and variant calling debugging tools')
def debug_tools():
  pass


@debug_tools.group('variant-call-analysis', short_help="Characterize TP, FN, FP and GT calls (and hence P/R) by variant size")
def variant_call_analysis():
  pass


@variant_call_analysis.command('process', short_help="Process EVCF, characterize calls by variant size and plot")
@click.argument('evcf', type=click.Path(exists=True))
@click.argument('out', type=click.Path())
@click.option('--region-label', help='Name of high confidence region if desired')
@click.option('--max-size', type=int, default=50, help='Maximum size of variant to consider')
@click.option('--title', help='Title for the plot')
@click.option('--fig-file', type=click.Path(), help='If supplied, plot will be saved here')
@click.option('--plot-bin-size', type=int, help='Bin size (bp)')
def vc_process(evcf, out, region_label, max_size, title, fig_file, plot_bin_size):
  import mitty.benchmarking.evcfbysize as ebs
  data = ebs.main(evcf_fname=evcf, out_csv_fname=out,
                  max_size=max_size, high_confidence_region=region_label)
  ebs.plot(data, fig_fname=fig_file, bin_size=plot_bin_size, title=title)


@variant_call_analysis.command('plot', short_help="Plot P/R from existing data file.")
@click.argument('datafile', type=click.Path(exists=True))
@click.argument('fig-file', type=click.Path())
@click.option('--title', help='Title for the plot')
@click.option('--plot-bin-size', type=int, help='Bin size (bp)')
@click.option('--plot-range', type=int, help='Range (bp) of indels to show')
def vc_process(datafile, fig_file, title, plot_bin_size, plot_range):
  import mitty.benchmarking.evcfbysize as ebs
  data = ebs.np.loadtxt(datafile, skiprows=1, delimiter=',', dtype=[('TP', int), ('FN', int), ('GT', int), ('FP', int)])
  ebs.plot(data, fig_fname=fig_file, bin_size=plot_bin_size, plot_range=plot_range, title=title)


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

  if fig_file is not None:
    vsd.plot(data, fig_fname=fig_file, bin_size=plot_bin_size, title=title)


@debug_tools.command('call-fate', short_help="Tracks fate of TP, FN, FP .. between two eval VCFs")
@click.argument('vcfa', type=click.Path(exists=True))
@click.argument('vcfb', type=click.Path(exists=True))
@click.argument('vcfout', type=click.File('w'))
@click.argument('summaryout', type=click.File('w'))
@click.option('--region-label', help='Name of high confidence region if desired')
def call_fate(vcfa, vcfb, vcfout, summaryout, region_label):
  """This tool tracks the fate of every variant call across the two supplied files and divides
  the calls up into the following 12 transitions

  \b
  Improvements
  ------------
  FN -> TP
  FN -> GT
  GT -> TP
  FP -> N  (FP calls removed)

  \b
  Status quo
  ----------
  TP -> TP
  FN -> FN
  GT -> GT
  FP -> FP

  \b
  Regressions
  -----------
  TP -> FN
  TP -> GT
  GT -> FN
  N  -> FP (New FP calls)

  """
  import mitty.benchmarking.callfate as cf
  cf.main(fname_a=vcfa, fname_b=vcfb, vcf_out=vcfout, summary_out=summaryout, high_confidence_region=region_label)


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

@debug_tools.command('bam-to-truth', short_help="from input bam with mapping quality threshold to produce truth fastq and its long qname file")
@click.argument('bam_fpath_in',type=click.Path(exists=True))
@click.argument('mq_threshold',type=int)
@click.argument('sample_name')
@click.argument('output_prefix')
def bam_to_truth(bam_fpath_in,mq_threshold,sample_name,output_prefix):
  """
  Given bam file and mapping quality threshold ,the tool outputs 2 fastqs with their longqnames 
  if both paired mate and the read have mq above.Also adds sample_name to output.
  For qname format specification, check mitty Readme.MD documentation.
  """
  import mitty.benchmarking.bam_to_truth as btt
  btt.bam_to_truth(bam_fpath_in,mq_threshold,sample_name,output_prefix)


@debug_tools.group('alignment-analysis', short_help='Plot various alignment metrics from BAM')
def alignment_analysis():
  """Computes a three dimensional histogram of alignment metrics from a BAM of simulated reads

  \b
  The dimensions are:
    [0] Xd - alignment error  -max_xd, ... 0, ... +max_xd, wrong_chrom, unmapped
                              (2 * max_xd + 3)
    [1] MQ - mapping quality  0, ... max_MQ
                              (max_MQ + 1)
    [2] vlen - length of variant carried by read
                              Ref, < -max_vlen , -max_vlen, ... 0, ... +max_vlen, > +max_vlen
                              ( 2 * max_vlen + 1 + 2 + 1)"""
  pass


@alignment_analysis.command('process', short_help='Compute alignment metrics from BAM and plot')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('long_qname_file', type=click.Path(exists=True))
@click.argument('out', type=click.Path())
@click.option('--max-d', type=int, default=200, help='Range of d_err to consider')
@click.option('--max-size', type=int, default=50, help='Maximum size of variant to consider')
@click.option('--fig-prefix', type=click.Path(), help='If supplied, a series of plots will be saved with this prefix')
@click.option('--plot-bin-size', default=1, type=int, help='Bin size')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--processes', default=2, help='How many processes to use for computation')
def alignment_debug_plot(bam, long_qname_file, out, max_d, max_size, fig_prefix, plot_bin_size, strict_scoring, processes):
  """Computes 3D matrix of alignment metrics (d_err, MQ, v_size) saves it to a numpy array file and produces a set
  of summary figures"""

  # For automated workflows we sometimes have real (not simulated) data. For such workflows we simply
  # plot a placeholder figure to include in reports by passing a dummy long qname with the magic word
  # 'deadbeef' in it
  if open(long_qname_file, 'r').read().strip() == 'deadbeef':
    if fig_prefix is not None:
      import mitty.benchmarking.plot.placeholderfigure as phf
      phf.placeholder_figure('NOT SIMULATED READS', fig_prefix=fig_prefix)
    exit(0)

  import mitty.benchmarking.xmv as xmv
  xmv_mat = xmv.main(bam, sidecar_fname=long_qname_file,
                     max_xd=max_d, max_vlen=max_size, strict_scoring=strict_scoring,
                     processes=processes)
  xmv.save(xmv_mat, out)

  if fig_prefix is not None:
    xmv.plot_figures(xmv_mat, fig_prefix=fig_prefix, plot_bin_size=plot_bin_size)


@alignment_analysis.command('plot', short_help='Plot alignment metrics from existing data file')
@click.argument('datafile', type=click.Path(exists=True))
@click.argument('fig-prefix', type=click.Path())
@click.option('--plot-bin-size', default=1, type=int, help='Bin size')
def alignment_debug_plot(datafile, fig_prefix, plot_bin_size):
  """Produces a set of alignment metric summary figures from existing data file"""
  import mitty.benchmarking.xmv as xmv
  xmv_mat = xmv.np.load(datafile)
  xmv.plot_figures(xmv_mat, fig_prefix=fig_prefix, plot_bin_size=plot_bin_size)


@debug_tools.command('subset-bam', short_help="Subset a BAM based on d_err and variant size")
@click.argument('bamin', type=click.Path(exists=True))
@click.argument('sidecar', type=click.Path(exists=True))
@click.argument('bamout', type=click.Path())
@click.option('--d-range', type=(int, int), default=(-200, 200))
@click.option('--reject-d-range', is_flag=True, help='Reject reads inside the range instead of outside')
@click.option('--v-range', type=(int, int), default=(-200, 200))
@click.option('--reject-v-range', is_flag=True, help='Reject reads inside the range instead of outside')
@click.option('--reject-reads-with-variants', is_flag=True, help='Reject any reads carrying variants')
@click.option('--reject-reference-reads', is_flag=True, help='Reject reads with no variants')
@click.option('--do-not-index', is_flag=True, help='Do not index BAM file')
@click.option('--strict-scoring', is_flag=True, help="Don't consider breakpoints when scoring alignment")
@click.option('--processes', default=2, help='How many processes to use for computation')
def subset_bam(bamin, sidecar, bamout,
               d_range, reject_d_range,
               v_range, reject_v_range,
               reject_reads_with_variants, reject_reference_reads,
               strict_scoring,
               do_not_index, processes):
  """Produce a subset of an input BAM based on d_err and variant size"""
  import mitty.benchmarking.subsetbam as sub
  assert d_range[0] <= d_range[1], 'd_range error ({})'.format(d_range)
  assert v_range[0] <= v_range[1], 'v_range error ({})'.format(v_range)
  assert not (reject_reads_with_variants and reject_reference_reads), 'Can not reject both variant and reference reads'
  sub.main(bam_fname=bamin, sidecar_fname=sidecar, out_fname=bamout,
           d_range=d_range, reject_d_range=reject_d_range,
           v_range=v_range, reject_v_range=reject_v_range,
           reject_reads_with_variants=reject_reads_with_variants,
           reject_reference_reads=reject_reference_reads,
           strict_scoring=strict_scoring, do_not_index=do_not_index, processes=processes)


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
    'illumina': mitty.simulation.illumina,
  }.get(model['model_class'])

  return read_module, model


@cli.group('utils', short_help='Miscellaneous utilities')
def utils():
  pass


@utils.command('retruncate-qname', short_help='Given a BAM of FASTQ')
@click.argument('mainfile-in', type=click.Path(exists=True))
@click.argument('sidecar-in', type=click.Path(exists=True))
@click.argument('mainfile-out', type=click.Path())
@click.argument('sidecar-out', type=click.Path())
@click.option('--truncate-to', type=int, default=240)
@click.option('--file-type', type=click.Choice(['BAM', 'FASTQ']), help='If supplied overrides autodetection')
def retruncate_qname(mainfile_in, sidecar_in, mainfile_out, sidecar_out, truncate_to, file_type):
  """Given a FASTQ (or BAM) and a long-qnames file trancate the qnames to
whatever length we want and push the too-long qnames to the overflow
file. Also convert any old style qnames (not ending in a *) to new style
ones with a proper termination character"""
  import mitty.empirical.qnametruncate as qt
  qt.main(mainfile_in, sidecar_in, mainfile_out, sidecar_out, truncate_to=truncate_to, file_type=file_type)


@utils.command('vcf-complexity', short_help='Annotate variants with complexity measures')
@click.argument('vcfin', type=click.Path(exists=True))
@click.argument('vcfout', type=click.Path())
@click.argument('ref', type=click.Path(exists=True))
@click.argument('bg', type=click.Path(exists=True))
@click.option('--window-size', type=int, default=100, help='Window size of SE and LC computation')
def vcf_complexity(vcfin, vcfout, ref, bg, window_size):
  """Annotate variants with complexity measures"""
  import mitty.benchmarking.complexity as cplx
  cplx.vcf_complexity(
    vcf_in_fname=vcfin, vcf_out_fname=vcfout,
    ref_fname=ref, bg_fname=bg, window_size=window_size)


@utils.command('gc-cov')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-b', '--block-len', type=int, default=10000, help='Block size for GC/cov computation')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def gc_cov(bam, fasta, pkl, block_len, threads):
  """Calculate GC content vs coverage from a BAM. Save in pickle file"""
  import mitty.empirical.gc as megc
  megc.process_bam_parallel(bam, fasta, pkl, block_len=block_len, threads=threads)


@utils.command('bq')
@click.argument('bam', type=click.Path(exists=True))
@click.argument('pkl')
@click.option('-t', '--threads', type=int, default=1, help='Threads to use')
def sample_bq(bam, pkl, threads):
  """BQ distribution from BAM"""
  import mitty.empirical.bq as bbq
  bbq.process_bam_parallel(bam, pkl, threads=threads)


@utils.command('filter-eval-vcf', short_help='Split out the FP and FN from an eval.vcf')
@click.argument('vcfin', type=click.Path(exists=True))
@click.argument('outprefix')
def filter_eval_vcf(vcfin, outprefix):
  """Subset VCF for given sample, apply BED file and filter out complex variants
   making it suitable to use for read generation"""
  import mitty.benchmarking.filterevalvcf as fev
  fev.extract_fp_fn(vcfin, outprefix)


@utils.command('qname-stats', short_help='Given a FASTQ + side-car show us qname distribution')
@click.argument('fastq', type=click.Path(exists=True))
@click.argument('sidecar', type=click.Path(exists=True))
@click.option('--max-expected-qname', type=int, default=500, help='qname bin upper bound')
def qname_stats(fastq, sidecar, max_expected_qname):
  """A simple routine to load in a FASTQ file and give us the distribution of
  qname lengths, because I was curious"""
  import mitty.empirical.qnamestats as qs
  qname_count = qs.main(fastq_fname=fastq, qname_overflow_fname=sidecar,
                        max_expected_qname_length=max_expected_qname)
  for n, ql in enumerate(qname_count):
    print('{}, {}'.format(n, ql))