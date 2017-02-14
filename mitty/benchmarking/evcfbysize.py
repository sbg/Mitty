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

from mitty.lib.evcfparse import parse_line


def main(evcf_fname, out_csv_fname, plot_fname=None, max_size=50, bin_size=5):
  data = np.zeros(shape=(2 * max_size + 1 + 2), dtype=[('TP', int), ('FN', int), ('GT', int), ('FP', int)])
  offset, max_n = max_size + 1, 2 * max_size + 2
  mode = 'rb' if evcf_fname.endswith('bcf') else 'r'
  vcf_in = pysam.VariantFile(evcf_fname, mode)
  for v in vcf_in:
    c, s = parse_line(v)
    # if c == 'skip':
    #   print('TP with diff rep')
    #   print(v)
    if c is None:
      print(v)
    elif c is not 'skip':
      data[c][min(max_n, max(0, offset + s))] += 1

  np.savetxt(out_csv_fname, data, fmt='%d', delimiter=', ', header='TP, FN, GT, FP')


def plot_panels(ax, y, color='k', linestyle='-'):
  """Interpret y (an array of length 2 * N + 3) as follows:
    The indexes represent variant sizes from

    DEL len > N
    DEL len from N to 1
    SNPs
    INS len from 1 to N
    INS len > N

  Plot all this on axis ax, but make five separate segments as written above

  :param ax:
  :param y:
  :param fmt:
  :return:
  """
  N = int((y.shape[0] - 3) / 2)
  parts = [
    ('>DEL', [-N - 1]),
    ('DEL', list(range(-N, 0))),
    ('SNP', [0]),
    ('INS', list(range(1, N + 2))),
    ('>INS', [N + 2])
  ]

  ax.plot(parts[0][1], y[0], 'o', color=color, linestyle=linestyle)
  ax.plot(parts[1][1], y[1:N + 1], color=color, linestyle=linestyle)
  ax.plot(parts[2][1], y[N + 1], 'o', color=color, linestyle=linestyle)
  ax.plot(parts[3][1], y[N + 2:2 * N + 3], color=color, linestyle=linestyle)
  ax.plot(parts[4][1], y[0], 'o', color=color, linestyle=linestyle)


# Ignoring bin size for now
def plot(data, fname, bin_size=5):
  fig = plt.figure(figsize=(6, 11))
  plt.subplots_adjust(bottom=0.05, top=0.99)

  ax1 = plt.subplot(411)
  plot_panels(ax1, data['TP'] / (data['TP'] + data['FN']))

  plt.savefig(fname)