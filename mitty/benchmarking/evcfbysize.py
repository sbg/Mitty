"""Code to characterize the TP, FN and FP (and hence P/R) by variant size."""
import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.stats as ss

from mitty.lib.evcfparse import parse_line


def main(evcf_fname, out_csv_fname, max_size=50, high_confidence_region=None):
  data = np.zeros(shape=(2 * max_size + 1 + 2), dtype=[('TP', int), ('FN', int), ('GT', int), ('FP', int)])
  offset, max_n = max_size + 1, 2 * max_size + 2
  mode = 'rb' if evcf_fname.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(evcf_fname, mode)
  for v in vcf_in:
    if high_confidence_region is not None and high_confidence_region not in v.info.get('Regions', []):
      continue
    c, s = parse_line(v)
    if c is None:
      print(v)
    elif c is not 'skip':
      data[c][min(max_n, max(0, offset + s))] += 1

  np.savetxt(out_csv_fname, data, fmt='%d', delimiter=', ', header='TP, FN, GT, FP')
  return data


def bootstrap(p, n, ci=0.05):
  """Given a fraction correct (p) and total count (n) (both arrays) and a confidence interval
  return us two arrays l, h indicating the lower and higher confidence limits

  :param p:
  :param n:
  :param ci:
  :return: recarray with fields p, l and h
  """
  l_ci, h_ci = ci/2, 1 - ci/2
  # print([(_p,) + tuple(ss.binom.ppf([l_ci, h_ci], _n, _p)/_n) for _n, _p in zip(n, p)])
  return np.array([(_p,) + tuple(ss.binom.ppf([l_ci, h_ci], _n, _p)/_n) for _n, _p in zip(n, p)],
                  dtype=[('y', float), ('l', float), ('h', float)])


def bin_data_by_variant_size(num, den=None, bin_size=None):
  """Interpret num and den (arrays of length 2 * N + 3) as follows:

  The indexes represent variant sizes from
    DEL len > N
    DEL len from N to 1
    SNPs
    INS len from 1 to N
    INS len > N

  if den is supplied num/den is the probability of success and den is the total count of observations.
      den is used to compute confidence intervals
  if den is None, num is considered a total count. No CI is computed

  Return a list of 5 dicts. Each dict contains information that helps us to position and label the data as well as
  the actual data itself.

  :param num: an array of length 2 * N + 3
  :param den: an array of length 2 * N + 3
  :param bin_size: If bin_size is None then return data with original granularity
  :return:
  """
  N_orig = int((num.shape[0] - 3) / 2)
  bin_size = bin_size or 1
  b_num = bin_array(num, b=bin_size)
  b_den = None if den is None else bin_array(den, b=bin_size)
  y = np.array(b_num, dtype=[('y', int)]) if b_den is None else bootstrap(b_num / (b_den + 1e-6), b_den)
  N = int((y.shape[0] - 3) / 2)
  return [
    {'x': [-N - 1], 'xticks': [-N - 1], 'xticklabels': ['DEL > {}'.format(N_orig)], 'y': y[0:1]},
    {'x': list(range(-N, 0)), 'xticks': [-(N + 1) / 2], 'xticklabels': ['DEL'], 'y': y[1:N + 1]},
    {'x': [0], 'xticks': [0], 'xticklabels': ['SNP'], 'y': y[N + 1:N + 2]},
    {'x': list(range(1, N + 1)), 'xticks': [(N + 1)/ 2], 'xticklabels': ['INS'], 'y': y[N + 2:2 * N + 2]},
    {'x': [N + 1], 'xticks': [N + 1], 'xticklabels': ['INS > {}'.format(N_orig)], 'y': y[2 * N + 2:]}
  ]


def bin_array(y, b=20):
  N = int((y.shape[0] - 3) / 2)
  q = int(np.ceil(N / b)) - 1 # number of bins
  # The only trick here is to consolidate the 1:N+1 and N+2:2*N+2 sections into bins
  return np.concatenate((
    y[0:1],
    [y[max(1, N - (i + 1) * b + 1):N - i * b + 1].sum() for i in range(q -1, -1, -1)],
    y[N + 1:N + 2],
    [y[N + 2 + i * b:min(2 * N + 1, N + 2 + (i + 1) * b)].sum() for i in range(q)],
    y[2 * N + 2:]))


def line_with_ci(ax, part, color='k', linestyle='-'):
  ax.plot(part['x'], part['y']['y'], marker=None if len(part['x']) > 1 else 'x', color=color, linestyle=linestyle)
  if 'l' in part['y'].dtype.names:
    if len(part['y']['l']) > 1:
      ax.fill_between(part['x'], part['y']['l'], part['y']['h'], color=color, alpha=0.5)
    else:
      ax.plot([part['x'][0], part['x'][0]], [part['y']['l'][0], part['y']['h'][0]], color=color, linestyle=linestyle, lw=3)


def plot_panels(ax, num, den=None, bin_size=None,
                color='k', linestyle='-', label='',
                yscale='linear', yticks=None, ylim=None, show_ax_label=False,
                ci=0.05):
  """Interpret y (an array of length 2 * N + 3) as follows:
    The indexes represent variant sizes from

    DEL len > N
    DEL len from N to 1
    SNPs
    INS len from 1 to N
    INS len > N

  :param ax:
  :param y: success rate
  :param cnt: number of counts in bin
  :param color:
  :param linestyle:
  :param label:
  :param yscale:
  :param yticks:
  :param ylim:
  :param show_ax_label:
  :return:
  """
  parts = bin_data_by_variant_size(num, den, bin_size)
  x1, x0 = parts[-1]['x'][0], parts[0]['x'][0]
  x_lim = [x0 - .1 *(x1 - x0), x1 + .1 *(x1 - x0)]
  xticks, xticklabels = [], []
  for part in parts:
    # print(part)
    line_with_ci(ax, part, color=color, linestyle=linestyle)
    xticks += part['xticks']
    xticklabels += part['xticklabels']

  # The fixed reference lines
  ax.axvline(x=parts[0]['x'], color='k', linestyle=':')
  ax.axvline(x=parts[2]['x'], color='k', linestyle=':')
  ax.axvline(x=parts[4]['x'], color='k', linestyle=':')

  ax.set_xticks(xticks)
  ax.set_xticklabels(xticklabels if show_ax_label else [])

  ax.set_yscale(yscale)
  if yticks is not None:
    ax.set_yticks(yticks)
  if ylim is not None:
    ax.set_ylim(ylim)
  ax.get_yaxis().get_major_formatter().labelOnlyBase = False
  ax.set_xlim(x_lim)

  return mpatches.Patch(color=color, linestyle=linestyle, label=label)  # patch_for_legend


# Ignoring bin size for now
def plot(data, fig_fname, bin_size=5):
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99, hspace=0.01)

  ax1 = plt.subplot(411)
  r_p = plot_panels(ax1,
                    num=data['TP'],
                    den=(data['TP'] + data['FN'] + data['GT']),
                    bin_size=bin_size,
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='b', label='recall')
  p_p = plot_panels(ax1,
                    num=data['TP'],
                    den=(data['TP'] + data['FP']),
                    bin_size=bin_size,
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='r', label='precision')
  plt.legend(handles=[r_p, p_p], loc='lower center', fontsize=9)

  ax2 = plt.subplot(412)
  gt_p = plot_panels(ax2,
                     num=data['GT'],
                     den=data['TP'],
                     bin_size=bin_size,
                     yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                     color='k', label='GT')
  plt.legend(handles=[gt_p], loc='upper left', fontsize=9)

  ax3 = plt.subplot(413)
  n_max = 10 ** np.ceil(np.log10(data['FP'].max()))
  n_min = 10 ** int(np.log10(data['FP'].min() + 1))
  fp_p = plot_panels(ax3,
                     num=data['FP'],
                     bin_size=bin_size,
                     color='k', label='FP',
                     yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                     yscale='log')
  plt.legend(handles=[fp_p], loc='upper right', fontsize=9)

  ax4 = plt.subplot(414)
  n_max = 10 ** np.ceil(np.log10((data['TP'] + data['FN'] + data['GT']).max()))
  n_min = 10 ** int(np.log10((data['TP'] + data['FN'] + data['GT']).min() + 1))
  tot_p = plot_panels(ax4,
                      num=data['TP'] + data['FN'] + data['GT'],
                      bin_size=bin_size,
                      yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                      color='k', label='Total', show_ax_label=True, yscale='log')
  plt.legend(handles=[tot_p], loc='upper right', fontsize=9)

  plt.savefig(fig_fname)