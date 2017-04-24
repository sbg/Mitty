import numpy as np
from matplotlib import patches as mpatches
from scipy import stats as ss


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
  return np.array([((_p,) + tuple(ss.binom.ppf([l_ci, h_ci], _n, _p)/_n)) if _p > 0 else (0, 0, 0) for _n, _p in zip(n, p)],
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
  bin_size = min(bin_size or 1, int(N_orig / 2.0))
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


def censor(y):
  """Use clip to turn zeros into tiny positive values"""
  return y.clip(1e-6, y.max())


def line_with_ci(ax, part, color='k', linestyle='-'):
  ax.plot(part['x'], censor(part['y']['y']), marker=None if len(part['x']) > 1 else 'x', color=color, linestyle=linestyle)
  if 'l' in part['y'].dtype.names:
    if len(part['y']['l']) > 1:
      ax.fill_between(part['x'], censor(part['y']['l']), censor(part['y']['h']), color=color, alpha=0.5)
    else:
      ax.plot([part['x'][0], part['x'][0]], [censor(part['y']['l'][0]), censor(part['y']['h'][0])], color=color, linestyle=linestyle, lw=3)


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