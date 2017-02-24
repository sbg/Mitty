"""Program to go through a VCF and compute the distribution of variant sizes"""
import numpy as np
import pysam

from mitty.benchmarking.plot.byvsize import plot_panels

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(vcf_fname, sample_name=None, max_size=50):
  data = np.zeros(shape=(2 * max_size + 1 + 2), dtype=int)
  offset, max_n = max_size + 1, 2 * max_size + 2
  mode = 'rb' if vcf_fname.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(vcf_fname, mode)
  if sample_name is not None:
    vcf_in.subset_samples([sample_name])
  for v in vcf_in:
    for alt in v.alts:
      sz = len(alt) - len(v.ref)
      data[min(max_n, max(0, offset + sz))] += 1
  return data


def plot(data, fig_fname, bin_size=5, title='Variant size distribution'):
  fig = plt.figure(figsize=(11, 6))
  plt.subplots_adjust(bottom=0.05, top=0.95, hspace=0.01)

  ax1 = plt.subplot(111)
  n_max = 10 ** np.ceil(np.log10(data.max()))
  n_min = 10 ** int(np.log10(data.min() + 1))
  r_p = plot_panels(ax1,
                    num=data,
                    bin_size=bin_size,
                    yticks=[n_min, n_max], ylim=[n_min / 2, n_max * 2],
                    yscale='log',
                    color='k',
                    label='Count',
                    show_ax_label=True)
  plt.suptitle(title)
  plt.savefig(fig_fname)