import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plot_read_counts(ax, counts_l, keys=None, colors=None, labels=None):
  """Plot a horizontal histogram of read counts.

  :param ax:
  :param counts_l: dictionary of counts.
                   If a list of dictionaries is passed multiple bars are plotted.
                   Dictionary is in the same format as returned e.g. by read_counts


  :param keys: List of keys. If supplied only these keys and in this order will be shown
  :param colors: if passed will control the color of the bars
  :param labels: if passed will supply a legend with these labels
  :return:
  """
  colors = colors or ['g', 'c', 'r', 'y', 'b']

  if not isinstance(counts_l, list):
    counts_l = [counts_l]
  dh = 0.8 / len(counts_l)
  all_keys = keys or set(k for count in counts_l for k in count.keys())
  for i, counts in enumerate(counts_l):
    for n, k in enumerate(all_keys):
      cnt = counts.get(k, 0)
      ax.barh(n + i * dh, cnt, color=colors[i], align='center', height=0.8 * dh)
      if cnt > 0: ax.text(cnt, n + i * dh, cnt, color='w', va='center', ha='right')

  y_ticks = list(range(len(all_keys)))
  ax.set_yticks(y_ticks)
  ax.set_yticklabels(all_keys)
  ax.set_ylim([y_ticks[0] - 0.5, y_ticks[-1] + 1])
  ax.set_ylabel('Read fate')

  ax.set_xscale('log')
  ax.set_xlabel('Read count')

  if labels is not None:
    l_p = [
      mpatches.Patch(color=colors[k], label=labels[k])
      for k in range(len(labels))
    ]
    plt.legend(handles=l_p, loc='best')


def plot_mean_MQ_vs_derr(ax, dmv_mat, fmt='y', ms=3):
  """Plot mean_MQ (y-axis) against d_err (x-axis) for given range of variant sizes

  :param dmv_mat:
  :param fig_prefix:
  :param plot_bin_size:
  :return:
  """
  mq_mat = dmv_mat.sum(axis=2)
  data_cnt = mq_mat.sum(axis=1)
  mq_vector = np.arange(mq_mat.shape[1])
  mean_mq = np.dot(mq_vector, mq_mat.T) / np.clip(data_cnt, 1, data_cnt.max() + 1).astype(float)

  max_d = int((dmv_mat.shape[0] - 3) / 2)
  x1 = max_d * 1.3
  xt = [-max_d, -int(max_d/2), 0, int(max_d/2), max_d]
  xlim = [-max_d * 1.1, x1 * 1.1]

  ax.plot(range(-max_d, max_d + 1), mean_mq[:2 * max_d + 1], fmt, ms=ms)
  ax.plot([x1], mean_mq[2 * max_d + 1:2 * max_d + 2], fmt, ms=ms)

  ax.set_xticks(xt + [x1])
  ax.set_xticklabels(xt + ['WC'])
  ax.set_xlim(xlim)

  ax.set_yticks([0, 20, 40, 60])
  ax.set_ylim([-5, 70])
  ax.axvline(x=0, color='k', linestyle=':')

  plt.xlabel(r'$d_{err}$')
  plt.ylabel('Mean MQ')


#TODO: remove joint plotting of two d_err thresholds - pass this as a paramter so we can compose the plots
def plot_perr_vs_MQ(ax, dmv_mat, yscale='linear'):
  mq = np.arange(dmv_mat.shape[1])
  mq_mat = dmv_mat.sum(axis=2)
  for fmt, d_margin, label in zip(['g*', 'ro'], [0, 50], [r'$d_{err} = 0$', r'$|d_{err}| \leq 50$']):
    p_err = compute_p_err(mq_mat, d_margin)
    ax.plot(mq, p_err, fmt, mfc='none', label=label)

  ax.plot(mq, 10 ** (-mq / 10), 'k:')
  plt.setp(ax, xlabel='MQ', ylabel=r'$p_{err}$', xlim=(-1, 65), yscale=yscale)
  # ax.legend(loc='upper right' if yscale=='linear' else 'lower left')


def compute_p_err(mq_mat, d_margin=0):
  d_max = int((mq_mat.shape[0] - 3)/2)
  c = mq_mat[d_max - d_margin:d_max + 1 + d_margin, :].sum(axis=0)
  t = mq_mat.sum(axis=0)
  return (t - c) / np.clip(t, 1, 1 + t.max())


def plot_read_count_vs_MQ(ax, dmv_mat):
  mq = np.arange(dmv_mat.shape[1])
  mq_mat = dmv_mat.sum(axis=2)
  cnt = mq_mat.sum(axis=0)
  ax.step(mq, cnt, 'k', where='mid')
  plt.setp(ax, xlabel='MQ', ylabel=r'Read count', xlim=(-1, 65), yscale='log')


def hengli_plot(ax, dmv_mat, v_bin_label, d_err=20, fmt='r.-', label='aligner',
                xlim=(1e-7, 1e-1), ylim=(0.95, 1.002)):
  """
    x-axis # wrong / # mapped
    y-axis # mapped / total reads

    Parametrized by MQ >= q where q ranges from 0 to 70

  :param ax:
  :param dmv_mat:
  :param v_bin_label:
  :param d_err:
  :return:
  """
  max_vlen = int((dmv_mat.shape[2] - 3) / 2)
  assert max_vlen >= 50

  v_bins = {
    'Ref': [-2, -1],
    'SNP': [max_vlen, max_vlen + 1],
    'INS <= 10': [max_vlen + 1, max_vlen + 1 + 10],
    '10 < INS <= 50': [max_vlen + 11, max_vlen + 11 + 50],
    'DEL <= 10': [max_vlen - 10, max_vlen],
    '10 < DEL <= 50': [max_vlen - 50, max_vlen - 10]
  }

  # max_mq = dmv_mat.shape[1] - 1
  max_mq = 60

  x = np.empty(shape=(max_mq,))
  y = np.empty(shape=(max_mq,))

  max_xd = int((dmv_mat.shape[0] - 3) / 2)

  vb = v_bins[v_bin_label]

  total = dmv_mat[:, :, vb[0]:vb[1]].sum()
  # print(v_bin_label, total)

  for m in range(max_mq):
    wrong = dmv_mat[0: max_xd - d_err, m:, vb[0]:vb[1]].sum() + \
            dmv_mat[max_xd + d_err:-1, m:, vb[0]:vb[1]].sum()
    mapped = dmv_mat[0:-1, m:, vb[0]:vb[1]].sum()

    x[m] = wrong / (mapped + 1e-6)
    y[m] = mapped / (total + 1e-6)

  ax.semilogx(x, y, fmt, label=label, lw=3, ms=3)
  ax.text(xlim[0] * 2, ylim[1] * 0.998, v_bin_label, ha='left', va='top',
          size=8)
  plt.grid(False)
  for p in ['top', 'right']:
    ax.spines[p].set(visible=False)

  ax.get_xaxis().tick_bottom()
  ax.get_yaxis().tick_left()
  ax.tick_params(which='both', direction='out')

  for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(9)

  plt.setp(ax, xlim=xlim, ylim=ylim)
