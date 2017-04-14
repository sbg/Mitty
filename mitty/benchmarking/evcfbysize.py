"""Code to characterize the TP, FN and FP (and hence P/R) by variant size."""
import matplotlib
import numpy as np
import pysam

from mitty.benchmarking.plot.byvsize import plot_panels

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mitty.lib.evcfparse import parse_line


def main(evcf_fname, out_csv_fname, max_size=50, high_confidence_region=None):
  data = np.zeros(shape=(2 * max_size + 1 + 2), dtype=[('TP', int), ('FN', int), ('GT', int), ('FP', int)])
  offset, max_n = max_size + 1, 2 * max_size + 2
  mode = 'rb' if evcf_fname.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(evcf_fname, mode)
  for v in vcf_in:
    if high_confidence_region is not None and high_confidence_region not in v.info.get('Regions', []):
      continue
    c, s, _ = parse_line(v)
    if c is None:
      print(v)
    elif c is not 'skip':
      data[c][min(max_n, max(0, offset + s))] += 1

  np.savetxt(out_csv_fname, data, fmt='%d', delimiter=', ', header='TP, FN, GT, FP')
  return data


# Ignoring bin size for now
def plot(data, fig_fname, bin_size=5, plot_range=None, title='P/R by variant size'):
  max_size = int((data.shape[0] - 3)/2)
  if plot_range is not None and plot_range < max_size:
    n0, n1 = max_size + 1 - (plot_range + 1), max_size + 1 + (plot_range + 1) + 1
    _data = data[n0:n1]
    for k in data.dtype.names:
      _data[0][k] = data[:n0][k].sum()
      _data[-1][k] = data[n1:][k].sum()
    data = _data

  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.95, hspace=0.01)

  ax1 = plt.subplot(411)
  r_p = plot_panels(ax1,
                    num=data['TP'] + data['GT'],
                    den=(data['TP'] + data['FN'] + data['GT']),
                    bin_size=bin_size,
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='b', label='recall')
  p_p = plot_panels(ax1,
                    num=data['TP'] + data['GT'],
                    den=(data['TP'] + data['FP'] + data['GT']),
                    bin_size=bin_size,
                    yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                    color='r', label='precision')
  plt.legend(handles=[r_p, p_p], loc='center left',fontsize=9)

  ax2 = plt.subplot(412)
  gt_p = plot_panels(ax2,
                     num=data['GT'],
                     den=data['TP'] + data['GT'],
                     bin_size=bin_size,
                     yticks=[0.0, 0.5, 1.0], ylim=[-0.05, 1.05],
                     color='k', label='GT')
  plt.legend(handles=[gt_p], loc='upper left', fontsize=9)

  ax3 = plt.subplot(413)
  n_max = 10 ** np.ceil(np.log10(data['FP'].max() + 1))
  n_min = 10 ** int(np.log10(data['FP'].min() + 1))
  fp_p = plot_panels(ax3,
                     num=data['FP'],
                     bin_size=bin_size,
                     color='k', label='FP',
                     yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                     yscale='log')
  plt.legend(handles=[fp_p], loc='upper right', fontsize=9)

  ax4 = plt.subplot(414)
  n_max = 10 ** np.ceil(np.log10((data['TP'] + data['FN'] + data['GT']).max() + 1))
  n_min = 10 ** int(np.log10((data['TP'] + data['FN'] + data['GT']).min() + 1))
  tot_p = plot_panels(ax4,
                      num=data['TP'] + data['FN'] + data['GT'],
                      bin_size=bin_size,
                      yticks=[n_min, n_max], ylim=[n_min/2, n_max * 2],
                      color='k', label='Total', show_ax_label=True, yscale='log')
  plt.legend(handles=[tot_p], loc='upper right', fontsize=9)
  plt.suptitle(title)

  plt.savefig(fig_fname)