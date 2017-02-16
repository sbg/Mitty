"""Code to characterize the TP, FN and FP (and hence P/R) by variant size.

mitty -v4 debug pr-by-size \
  evcf.in.vcf.gz \
  out.csv \
  --max-size 1000 \
  --bin-size 20 \
  --plot pr.size.pdf


"""
import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.stats as ss

from mitty.lib.evcfparse import parse_line


def main(evcf_fname, out_csv_fname, plot_fname=None, max_size=50, bin_size=5, high_confidence_region=None):
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
  :return: recarray with fields l and h
  """
  l_ci, h_ci = ci/2, 1 - ci/2
  return np.array([tuple(ss.binom.ppf([l_ci, h_ci], _n, _p)/_n) for _n, _p in zip(n, p)], dtype=[('l', float), ('h', float)])


def line_with_ci(ax, x, y, cnt=None, fmt='', color='k', linestyle='-'):
  ax.plot(x, y, fmt, color=color, linestyle=linestyle)
  if cnt is not None:
    ci_lh = bootstrap(y, cnt)
    if len(x) > 1:
      ax.fill_between(x, ci_lh['l'], ci_lh['h'], color=color, alpha=0.5)
    else:
      ax.plot([x[0], x[0]], [ci_lh['l'][0], ci_lh['h'][0]], color=color, linestyle=linestyle, lw=3)


def plot_panels(ax, y, cnt=None, color='k', linestyle='-', label='',
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
  q = not cnt is None

  N = int((y.shape[0] - 3) / 2)
  parts = [
    [-N - 1],
    list(range(-N, 0)),
    [0],
    list(range(1, N + 1)),
    [N + 1]
  ]

  line_with_ci(ax, parts[0], y[0:1], cnt=cnt[0:1] if q else None, fmt='x', color=color, linestyle=linestyle)
  line_with_ci(ax, parts[1], y[1:N + 1], cnt=cnt[1:N + 1] if q else None, color=color, linestyle=linestyle)
  line_with_ci(ax, parts[2], y[N + 1:N + 2], cnt=cnt[N + 1:N + 2] if q else None, fmt='x', color=color, linestyle=linestyle)
  line_with_ci(ax, parts[3], y[N + 2:2 * N + 2], cnt=cnt[N + 2:2 * N + 2] if q else None, color=color, linestyle=linestyle)
  line_with_ci(ax, parts[4], y[2 * N + 2:], cnt=cnt[2 * N + 2:] if q else None, fmt='x', color=color, linestyle=linestyle)

  # The fixed reference lines
  ax.axvline(x=parts[0], color='k', linestyle=':')
  ax.axvline(x=parts[2], color='k', linestyle=':')
  ax.axvline(x=parts[4], color='k', linestyle=':')

  # The special xtick labels
  if show_ax_label:
    ax.set_xticks([-N - 1, -N/2, 0, N/2, N + 1])
    ax.set_xticklabels(['DEL > {}'.format(N), 'DEL', 'SNP', 'INS', 'INS > {}'.format(N)])
  else:
    ax.set_xticklabels([])
  # plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)

  ax.set_yscale(yscale)
  if yticks is not None:
    ax.set_yticks(yticks)
  if ylim is not None:
    ax.set_ylim(ylim)
  ax.get_yaxis().get_major_formatter().labelOnlyBase = False

  return mpatches.Patch(color=color, linestyle=linestyle, label=label)  # patch_for_legend


# Ignoring bin size for now
def plot(data, fname, bin_size=5):
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99, hspace=0.01)

  ax1 = plt.subplot(411)
  r_p = plot_panels(ax1,
                    y = data['TP'] / (data['TP'] + data['FN'] + data['GT'] + 1e-6),
                    cnt = (data['TP'] + data['FN'] + data['GT']),
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='b', label='recall')
  p_p = plot_panels(ax1,
                    y = data['TP'] / (data['TP'] + data['FP'] + 1e-6),
                    cnt = (data['TP'] + data['FP']),
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='r', label='precision')
  plt.legend(handles=[r_p, p_p], loc='lower center', fontsize=9)

  ax2 = plt.subplot(412)
  gt_p = plot_panels(ax2,
                     y = data['GT'] / (data['TP'] + 1e-6),
                     cnt = data['TP'],
                     yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                     color='k', label='GT')
  plt.legend(handles=[gt_p], loc='upper left', fontsize=9)

  ax3 = plt.subplot(413)
  n_max = 10 ** np.ceil(np.log10(data['FP'].max()))
  n_min = 10 ** int(np.log10(data['FP'].min() + 1))
  fp_p = plot_panels(ax3,
                     y = data['FP'],
                     color='k', label='FP',
                     yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                     yscale='log')
  plt.legend(handles=[fp_p], loc='upper right', fontsize=9)

  ax4 = plt.subplot(414)
  n_max = 10 ** np.ceil(np.log10((data['TP'] + data['FN'] + data['GT']).max()))
  n_min = 10 ** int(np.log10((data['TP'] + data['FN'] + data['GT']).min() + 1))
  tot_p = plot_panels(ax4, data['TP'] + data['FN'] + data['GT'],
                      yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                      color='k', label='Total', show_ax_label=True, yscale='log')
  plt.legend(handles=[tot_p], loc='upper right', fontsize=9)

  plt.savefig(fname)