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