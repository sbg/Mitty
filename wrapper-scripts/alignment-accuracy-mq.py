"""
Script that takes in an aaf file and extracts mapping quality and alignment
accuracy.


Commandline is of the form

  python3 alignment-accuracy-mq.py  <AAF> <OUTNAME> <NCPUS>


Data is in the form of a csv that has the following columns

variant category, MQ, d_err_bin, read_count
"""
import sys

import mitty.analysis as maly

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


def count_alignment_category(riter):
  v_cats = ['Ref', 'Multi', 'SNP', 'DEL <= 10', 'DEL 11-50', 'INS <= 10', 'INS 11-50', 'DEL > 50', 'INS > 50']
  d_cats = ['d = 0', '0 < d <= 50', 'd > 50', 'wrong chrom', 'unmapped']

  un_cat = 0
  cnt = [[[0 for _ in v_cats] for _ in range(61)] for _ in d_cats]

  for r in riter:
    d_err, mq, v_list = r

    if len(v_list) == 0:
      i = 0  # ['Ref']
    elif len(v_list) > 1:
      i = 1  # ['Multi']
    elif v_list[0] == 0:
      i = 2  # ['SNP']
    elif -10 <= v_list[0] < 0:
      i = 3  # ['DEL <= 10']
    elif -50 <= v_list[0] < -10:
      i = 4  # ['DEL 11-50']
    elif 0 < v_list[0] <= 10:
      i = 5  # ['INS <= 10']
    elif 10 < v_list[0] <= 50:
      i = 6  # ['INS 11-50']
    elif v_list[0] < -50:
      i = 7  # ['DEL > 50']
    elif 50 < v_list[0]:
      i = 8  # ['INS > 50']
    else:
      un_cat += 1  # This should be zero, or we have a bug
      continue

    j = min(mq, 60)

    if d_err == 0:
      k = 0  # ['0']
    elif abs(d_err) <= 50:
      k = 1  # ['<= 50']
    elif 50 < abs(d_err) < 201:
      k = 2  # ['> 50']
    elif d_err == 201:
      k = 3  # ['wrong chrom']
    else:
      k = 4  # ['unmapped']

    cnt[k][j][i] += 1

  logger.debug('{} reads out of variant range'.format(un_cat))

  counter_dict = {}
  for i, v_cat in enumerate(v_cats):
    for j in range(61):
      for k, d_cat in enumerate(d_cats):
        counter_dict[','.join([v_cat, str(j), d_cat])] = cnt[k][j][i]

  return counter_dict


def combine_counts(counter_dict_l):
  counter_dict = None
  for cd in counter_dict_l:
    if counter_dict is None:
      counter_dict = {k: v for k, v in cd.items()}
    else:
      for k, v in cd.items():
        counter_dict[k] += v
  return counter_dict


if len(sys.argv) != 4:
  print('Usage: ')
  print('python3 {} <AAF> <OUTNAME> <NCPUS>'.format(sys.argv[0]))
  exit(1)

aaf_fname = sys.argv[1]
out_fname = sys.argv[2]
ncpus = int(sys.argv[3])

n1 = count_alignment_category
pipeline = [n1]

counts = combine_counts(maly.aaftoolz.scatter_aaf(pipeline, aaf_fname=aaf_fname,
                                         ncpus=ncpus))

with open(out_fname, 'w') as fout:
  fout.write(','.join(['category', 'MQ', 'derr', 'count']) + '\n')
  for k, v in counts.items():
    fout.write('{},{}\n'.format(k, v))