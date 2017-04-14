"""Takes the variants from two eVCF files and tracks the fate of each variant's category
between the two files."""
import time

import numpy as np
import pysam

from mitty.lib.evcfparse import parse_line

import logging
logger = logging.getLogger(__name__)


def create_partition_list():
  partitions = {
    '{}-{}'.format(a, b): {'count': 0, 'fp': None}
    for a in ['TP', 'FN', 'GT']
    for b in ['TP', 'FN', 'GT']
  }

  partitions.update(
    {'{}-{}'.format(a, b): {'count': 0, 'fp': None}
     for a in ['FP', 'N']
     for b in ['FP', 'N']
     if not (a == b == 'N')
     }
  )

  return partitions


def key_v(v, alt):
  return '{}.{}.{}.{}'.format(v.chrom, v.pos, v.ref, alt)


def process_v(v):
  cat, sz, alt = parse_line(v) if v is not None else (None, None, None)
  k = key_v(v, alt) if sz is not None else None
  return k, cat


def update(v, var_k, var_cat, file_no, partitions, working_dict):
  """
  :param v: the variant
  :param var_k:
  :param var_cat:
  :param file_no: 0 or 1 depending on which file the variant came from
  :param partitions:
  :param working_dict:
  :return:
  """
  val = working_dict.get(var_k, None)
  if val is not None:  # You complete me
    if val[1 - file_no] is None:
      raise ValueError('Repeated identical entry: {}'.format(var_k))
    val[file_no] = {'cat': var_cat, 'v': v}
    cat_key = '{}-{}'.format(val[0]['cat'], val[1]['cat'])
    partitions[cat_key]['count'] += 1
    partitions[cat_key]['vcf'].write(val[0]['v'])
    del working_dict[var_k]
  elif var_k is not None:
    working_dict[var_k] = [None, None]
    working_dict[var_k][file_no] = {'cat': var_cat, 'v': v}


def main(fname_a, fname_b, out_prefix, high_confidence_region=None):
  partitions = create_partition_list()
  vcf_a = pysam.VariantFile(fname_a, mode='rb' if fname_a.endswith('bcf') else 'r')
  vcf_b = pysam.VariantFile(fname_b, mode='rb' if fname_a.endswith('bcf') else 'r')
  for k, v in partitions.items():
    v['vcf'] = pysam.VariantFile(out_prefix + '-' + k + '.vcf', mode='w', header=vcf_a.header)

  logger.debug('Reading variants ...')
  t0 = time.time()
  ctr = 0
  working_dict = {}
  while 1:
    va, vb = next(vcf_a, None), next(vcf_b, None)
    if va is None and vb is None: break

    if high_confidence_region is not None:
      if high_confidence_region not in va.info.get('Regions', []):
        va = None
      if high_confidence_region not in vb.info.get('Regions', []):
        vb = None

    var_k_a, var_cat_a = process_v(va)
    var_k_b, var_cat_b = process_v(vb)

    update(va, var_k_a, var_cat_a, 0, partitions, working_dict)
    update(vb, var_k_b, var_cat_b, 1, partitions, working_dict)

    if var_k_a is not None: ctr += 1
    if var_k_b is not None: ctr += 1

    if 0 <= ctr % 100000 <= 1:
      t1 = time.time()
      logger.debug('Read {} variants ({:0.2f} v/s). Working dict size is {}'.format(ctr, ctr/(t1 - t0), len(working_dict)))

  t1 = time.time()
  logger.debug(
    'Finished reading {} variants ({:0.2f} v/s). Residual working dict size is {}'.format(ctr, ctr / (t1 - t0), len(working_dict)))

  # At this point, the working_dict better only contain N -> FP and FP -> N
  for k, v in working_dict.items():
    v_to_use, cat_k = (v[1], 'N-FP') if v[0] is None else (v[0], 'FP-N')
    assert v_to_use['cat'] == 'FP', '{}: Not a FP. There is a bug in the code'.format(k)
    partitions[cat_k]['count'] += 1
    partitions[cat_k]['vcf'].write(v_to_use['v'])

  summary = summary_table_text(partitions)
  logger.debug(summary)

  with open(out_prefix + '-summary.txt', 'w') as fp:
    fp.write(summary)


def summary_table_text(partitions):
  def line(k):
    return ['{}: {}'.format(k.replace('-', ' -> '), partitions[k]['count'])]

  res = ['', 'Improvements', '--------------']
  for key in ['FN-TP', 'FN-GT', 'GT-TP', 'FP-N']:
    res += line(key)

  res += ['', 'Unchanged', '--------------']
  for key in ['TP-TP', 'FN-FN', 'GT-GT', 'FP-FP']:
    res += line(key)

  res += ['', 'Regression', '--------------']
  for key in ['TP-FN', 'TP-GT', 'GT-FN', 'N-FP']:
    res += line(key)

  return '\n'.join(res)