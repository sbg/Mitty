"""Example read model

TODO: Figure out the logistics of having diffrent models and loading GC-bias, sequencing error etc.

"""
import numpy as np

SEED_MAX = (1 << 32) - 1  # Used for seeding rng


def read_model_params(model, diploid_coverage=30.0):
  """

  :param model
  :param diploid_coverage: Coverage we expect from a diploid genome
  :return: dict of params,
           passes -> how many times we should call generate reads for each region
  """
  rlen = model['mean_rlen']
  # coverage = read count * read len / genome len
  # coverage = P * read len P = expected number of reads per base
  # P = p * passes    p = probability of read per pass
  # One template for two reads
  # p = P / (2 * passes)
  # p = coverage / (2 * rlen * passes)
  p = 1.0
  passes = 1
  while p > 0.1:
    passes *= 2
    p = 0.5 * diploid_coverage / (2 * rlen * passes)

  return {
    'diploid_coverage': diploid_coverage,
    'p': p,
    'passes': passes,
    'rlen': rlen,
    'cum_tlen': model['cum_tlen'],
    'cum_bq_mat': model['cum_bq_mat'],
    'unpaired': model.get('unpaired', False)
  }


def generate_reads(model, p_min, p_max, seed=7):
  """

  :param model: dict as returned by read_model_params
  :param p_min: Start of reads
  :param p_max: End of reads
  :param region: ('1', 10001, 103863906) Like the BED file
  :param seed:
  :return:
  """
  if not (0 <= seed <= SEED_MAX):
    raise ValueError('Seed value {} is out of range 0 - {}'.format(seed, SEED_MAX))

  seed_rng = np.random.RandomState(seed)
  tloc_rng, tlen_rng, shuffle_rng, file_order_rng = \
    [np.random.RandomState(s) for s in seed_rng.randint(SEED_MAX, size=4)]

  return _reads_for_template_in_region(
        model, file_order_rng,
        *_templates_for_region(
          p_min, p_max, model, tloc_rng, tlen_rng, shuffle_rng))


def _templates_for_region(p_min, p_max, model, tloc_rng, tlen_rng, shuffle_rng):
  # TODO: GC-bias goes here
  p, rlen = model['p'], model['rlen']
  est_block_size = int((p_max - p_min) * p * 1.2)
  ts = tloc_rng.geometric(p=p, size=est_block_size).cumsum() + p_min + 1
  shuffle_rng.shuffle(ts)
  tl = np.searchsorted(model['cum_tlen'], tlen_rng.rand(ts.shape[0]))
  tl.clip(rlen, out=tl)
  te = ts + tl
  idx = (te < p_max)
  return ts[idx], te[idx]


def _reads_for_template_in_region(model, file_order_rng, ts, te):
  """

  :param rlen: Fixed - Illumina reads
  :param file_order_rng: Decides which one of the pair (F or R) comes first in the file (or goes in file1 of the pair)
  :param ts: Start of template
  :param te: End of template
  :return:
  """
  rlen = model['rlen']

  r0l = np.full(ts.size, rlen, dtype=np.uint32)
  r1l = np.full(ts.size, rlen, dtype=np.uint32)

  r0fo = file_order_rng.randint(2, size=ts.size, dtype='i1')
  r1fo = 1 - r0fo

  r0p = ts
  r1p = te - rlen

  if model.get('unpaired', False):
    return [
      {
        'file_order': np.zeros((2 * ts.size,), dtype='i1'),
        'strand': np.concatenate((np.full(ts.size, 0, dtype=np.uint8), np.full(ts.size, 1, dtype=np.uint8))),
        'pos': np.concatenate((r0p, r1p)),
        'len': np.concatenate((r0l, r1l))
      }
    ]
  else:
    return [
      {
        'file_order': r0fo,
        'strand': np.full(ts.size, 0, dtype=np.uint8),
        'pos': r0p,
        'len': r0l
      },
      {
        'file_order': r1fo,
        'strand': np.full(ts.size, 1, dtype=np.uint8),
        'pos': r1p,
        'len': r1l
      }
    ]


def corrupt_template(model, template, corrupt_rng):
  """Given a pair of read sequences, return a corrupted read pair

  :param model:
  :param template:
    (
      index,  - same as the original index from the perfect reads
      sample_name,
      chrom,
      copy,
      (
        (strand, pos, cigar, (v1,v2,...), MD, seq, qual)
        ...  [repeated as for as many reads in this template]
      )
    )

  :param corrupt_rng:
  :return:
    (
      index,  - same as the original index from the perfect reads
      sample_name,
      chrom,
      copy,
      (
        (strand, pos, cigar, (v1,v2,...), MD, seq, qual)
        ...  [repeated as for as many reads in this template]
      )
    )
  """
  bq_mat = model['cum_bq_mat']
  # index, sample_name, chrom, cpy = template[:4]
  return template[:4] + (tuple(
    ((strand, pos, cigar, v_list) + corrupt_single_read(seq, bq_mat[mate, :, :], corrupt_rng)
     for (strand, pos, cigar, v_list, _, seq, _), mate in zip(template[4], [0, 1]))
  ),)


base_rot = {
  'A': 'CTG',
  'C': 'ATG',
  'T': 'ACG',
  'G': 'ACT'
}  # We use this for creating base call errors
phred_p = 10 ** (-np.arange(100)/10)


def corrupt_single_read(seq, bq_mat, corrupt_rng):
  """Given a sequence, BQ model and rng, add errors to the sequence and return us a BQ quality array

  :param seq:
  :param bq_mat:
  :param corrupt_rng:
  :return: MD, seq, qual
  """
  rlen = len(seq)
  corrupt_seq = list(seq)
  bq_seq = [0] * rlen  # Just an initializer
  bq_rnd = corrupt_rng.rand(rlen)
  base_call_rnd = corrupt_rng.rand(rlen)
  base_rnd = corrupt_rng.randint(0, 3, size=rlen)

  md, ctr = '', 0
  for n in range(rlen):
    bq_seq[n] = bq = min(np.searchsorted(bq_mat[n, :], bq_rnd[n]), 93)
    # This min should only ever be invoked if our BQ model matrix is smaller than the read length
    # Since we extract read length and BQ model from the same BAM this should not happen
    if base_call_rnd[n] < phred_p[bq]:
      corrupt_seq[n] = base_rot.get(seq[n], 'NNN')[base_rnd[n]]
      md += str(ctr) + seq[n]
      ctr = 0
    else:
      ctr += 1

  return md, ''.join(corrupt_seq), ''.join([chr(b + 33) for b in bq_seq])


def describe_model(model_name, model, figfile):
  """Plot a few panels describing what the model looks like

  :param model_name:
  :param model:
  :param figfile:
  :return:
  """
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  fig = plt.figure(figsize=(6, 15))
  plt.subplots_adjust(bottom=0.05, top=0.92, hspace=0.4)
  fig.suptitle('MODEL: {}'.format(model_name), ha='left')

  ax = plt.subplot(4,1,1)
  plot_template_length_distribution(model, ax, plt)

  ax = plt.subplot(4,1,2)
  plot_BQ_heatmap(model, 0, ax, plt)

  ax = plt.subplot(4,1,3)
  plot_BQ_heatmap(model, 1, ax, plt)

  fig.text(0.1, 0.25, 'Model description:\n\n' + model['model_description'], va='top', wrap='True')

  plt.savefig(figfile)


def plot_template_length_distribution(model, ax, plt):
  ax.plot(model['tlen']/model['tlen'].max())
  ax.text(len(model['tlen']), plt.getp(ax, 'ylim')[1], 'n={}'.format(model['tlen'].sum()), ha='right', va='top')
  plt.setp(ax, xlabel='Template length (bp)', ylim=(0, 1.2), yticks=[], title='Template length distribution')


def plot_BQ_heatmap(model, mate, ax, plt):
  this_bq_mat = model['bq_mat'][mate, :model['max_rlen'], :65]
  plt.imshow(this_bq_mat.T,
             extent=(1, model['max_rlen'], 1, 65),
             aspect='auto',
             origin='lower', cmap=plt.cm.gray_r)
  # bq_mean = np.dot(this_bq_mat, np.arange(65)) / this_bq_mat.sum(axis=1)
  p_err = (np.dot(this_bq_mat, 10 ** (-np.arange(65)/10.0)) / this_bq_mat.sum(axis=1))
  bq_mean = -10 * np.log10(p_err)
  p_err_mean = p_err.mean()

  plt.plot(range(1, model['max_rlen'] + 1), bq_mean, 'y', lw=2)
  plt.text(1, 20, 'Mean error {:2.3}%'.format(p_err_mean * 100))
  plt.setp(ax, xlabel='Position on read (bp)', ylabel='BQ value', title='Base quality: Mate {}'.format(mate + 1))

