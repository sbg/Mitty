"""
Simple function to plot a small figure with a message we pass it
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def placeholder_figure(msg, fig_prefix):
  fig = plt.figure(figsize=(4, 4))
  plt.text(0.5, 0.5, msg, ha='center', va='center')
  ax = plt.gca()
  plt.setp(ax, xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[])
  plt.savefig(fig_prefix + '_placeholder.png')