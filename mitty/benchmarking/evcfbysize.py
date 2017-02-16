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
    # if c == 'skip':
    #   print('TP with diff rep')
    #   print(v)
    if c is None:
      print(v)
    elif c is not 'skip':
      data[c][min(max_n, max(0, offset + s))] += 1

  np.savetxt(out_csv_fname, data, fmt='%d', delimiter=', ', header='TP, FN, GT, FP')


def plot_panels(ax, y, color='k', linestyle='-', label='',
                yscale='linear', yticks=None, ylim=None, show_ax_label=False):
  """Interpret y (an array of length 2 * N + 3) as follows:
    The indexes represent variant sizes from

    DEL len > N
    DEL len from N to 1
    SNPs
    INS len from 1 to N
    INS len > N

  :param ax:
  :param y:
  :param color:
  :param linestyle:
  :param label:
  :param yscale:
  :param yticks:
  :param ylim:
  :param show_ax_label:
  :return:
  """
  N = int((y.shape[0] - 3) / 2)
  parts = [
    [-N - 1],
    list(range(-N, 0)),
    [0],
    list(range(1, N + 2)),
    [N + 2]
  ]

  ax.plot(parts[0], y[0], 'o', color=color, linestyle=linestyle)
  ax.axvline(x=parts[0], color='k', linestyle=':')
  ax.plot(parts[1], y[1:N + 1], color=color, linestyle=linestyle)
  ax.plot(parts[2], y[N + 1], 'o', color=color, linestyle=linestyle)
  ax.axvline(x=parts[2], color='k', linestyle=':')
  ax.plot(parts[3], y[N + 2:2 * N + 3], color=color, linestyle=linestyle)
  ax.plot(parts[4], y[0], 'o', color=color, linestyle=linestyle)
  ax.axvline(x=parts[4], color='k', linestyle=':')

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
  r_p = plot_panels(ax1, data['TP'] / (data['TP'] + data['FN'] + 1e-6),
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='b', label='recall')
  p_p = plot_panels(ax1, data['TP'] / (data['TP'] + data['FP'] + 1e-6),
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='r', label='precision')
  plt.legend(handles=[r_p, p_p], loc='lower center', fontsize=9)

  ax2 = plt.subplot(412)
  gt_p = plot_panels(ax2, data['GT'] / (data['TP'] + 1e-6),
                     yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                     color='k', label='GT')
  plt.legend(handles=[gt_p], loc='upper left', fontsize=9)

  ax3 = plt.subplot(413)
  n_max = 10 ** np.ceil(np.log10(data['FP'].max()))
  n_min = 10 ** int(np.log10(data['FP'].min() + 1))
  fp_p = plot_panels(ax3, data['FP'], color='k', label='FP',
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