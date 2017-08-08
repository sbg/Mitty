import numpy as np
import matplotlib.pyplot as plt


def plot_read_counts(ax, counts_l, keys=None, colors=None):
  """Plot a horizontal histogram of read counts.

  :param ax:
  :param counts_l: list of count dictionaries. All counts will be plotted on same axes
  :param keys: List of keys. If supplied only these keys and in this order will be shown
  :return:
  """
  colors = colors or ['g', 'c', 'r', 'y', 'b']

  if not isinstance(counts_l, list):
    counts_l = [counts_l]
  dh = 0.8 / len(counts_l)
  all_keys = keys or set(k for count in counts_l for k in count.keys())
  for i, counts in enumerate(counts_l):
    for n, k in enumerate(all_keys):
      ax.barh(n + i * dh, counts[k], color=colors[i], align='center', height= 0.8 * dh)
      ax.text(counts[k], n + i * dh, counts[k], color='w', va='center', ha='right')

  y_ticks = list(range(len(all_keys)))
  ax.set_yticks(y_ticks)
  ax.set_yticklabels(all_keys)
  ax.set_ylim([y_ticks[0] - 0.5, y_ticks[-1] + 1])
  ax.set_ylabel('Read fate')

  ax.set_xscale('log')
  ax.set_xlabel('Read count')


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
