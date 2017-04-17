"""Takes the variants from two eVCF files and tracks the fate of each variant's category
between the two files."""
import time

import numpy as np
import pysam

from mitty.lib.evcfparse import parse_line

import logging
logger = logging.getLogger(__name__)


def gt_str(s):
  return '{}{}{}'.format(s['GT'][0], '|' if s.phased else '/', s['GT'][1])


def partition_list():
  return [
    'FN->TP',
    'FN->GT',
    'GT->TP',
    'FP->N',
    'TP->TP',
    'FN->FN',
    'GT->GT',
    'FP->FP',
    'TP->FN',
    'TP->GT',
    'GT->FN',
    'N->FP'
  ]


def write_output_header(out_fp, vcf_fp):
  """Write the header for the output file

  :param out_fp:
  :param vcf_fp:
  :return:
  """
  hdr = vcf_fp.header.copy()
  for k in hdr.info.keys():
    if k != 'Regions':
      hdr.info.remove_header(k)
  for k in hdr.formats.keys():
    if k != 'GT':
      hdr.formats.remove_header(k)

  hdr_lines = str(hdr).splitlines(keepends=True)
  column_line = '\t'.join(hdr_lines[-1].split('\t')[:-2] + partition_list())  # Add the 12 samples

  hdr_lines[-1] = '##INFO=<ID=Vtype,Number=1,Type=String,Description="High-level variant type (SNP|INDEL).">\n'
  hdr_lines += [column_line]

  out_fp.writelines(hdr_lines)


def partition_count_dict():
  return {
    k: {vtype: 0 for vtype in ['SNP', 'INDEL']}
    for k in partition_list()
  }


def variant_hash(v, alt):
  return '{}.{}.{}.{}'.format(v.chrom, v.pos, v.ref, alt)


def process_v(v):
  cat, sz, alt = parse_line(v) if v is not None else (None, None, None)
  k = variant_hash(v, alt) if sz is not None else None
  return k, cat, sz


def construct_info_field(v, var_type):
  info_str = 'Vtype=' + var_type
  info = v.info.get('Regions', None)
  info_str += ';Regions={}'.format(','.join(info)) if info is not None else ''
  return info_str


def write_out_v(out_fp, v, var_type, var_change_cat):
  row = [str(v.chrom), str(v.pos), '.', v.ref, ','.join(v.alts), '.',
         'LowQual' if v.filter.get('LowQual', None) else '.', construct_info_field(v, var_type), 'GT']
  gt = gt_str(v.samples['QUERY']) \
    if var_change_cat[:2] == 'FP' or var_change_cat[-2:] == 'FP' \
    else gt_str(v.samples['TRUTH'])
  row += ['0/0' if k != var_change_cat else gt for k in partition_list()]
  out_fp.write('\n' + '\t'.join(row))


def update(out_fp, v, var_k, var_cat, var_type, file_no, partitions, working_dict):
  """
  :param out_fp:
  :param v: the variant
  :param var_k:
  :param var_cat:
  :param var_type: 'SNP' or 'INDEL'
  :param file_no: 0 or 1 depending on which file the variant came from
  :param partitions:
  :param working_dict:
  :return:
  """
  val = working_dict.get(var_k, None)
  if val is not None:  # You complete me
    if val[1 - file_no] is None:
      raise ValueError('Repeated identical entry: {}'.format(var_k))
    val[file_no] = {'cat': var_cat, 'type': var_type, 'v': v}
    cat_key = '{}->{}'.format(val[0]['cat'], val[1]['cat'])
    partitions[cat_key][var_type] += 1
    write_out_v(out_fp, val[0]['v'], var_type, cat_key)  # We use val[0] because these are guaranteed not to be N->FP
    del working_dict[var_k]
  elif var_k is not None:
    working_dict[var_k] = [None, None]
    working_dict[var_k][file_no] = {'cat': var_cat, 'type': var_type, 'v': v}


def main(fname_a, fname_b, vcf_out, summary_out, high_confidence_region=None):
  partitions = partition_count_dict()

  vcf_in = [
    pysam.VariantFile(fname_a, mode='rb' if fname_a.endswith('bcf') else 'r'),
    pysam.VariantFile(fname_b, mode='rb' if fname_a.endswith('bcf') else 'r')
  ]
  write_output_header(vcf_out, vcf_in[0])

  logger.debug('Reading variants ...')
  t0 = time.time()
  ctr, null_ctr = 0, 0
  working_dict = {}
  while null_ctr < 2:
    null_ctr = 0
    for file_no in [0, 1]:
      v = next(vcf_in[file_no], None)
      if v is None: null_ctr += 1
      if high_confidence_region is not None:
        if high_confidence_region not in v.info.get('Regions', []):
          v = None

      var_k, var_cat, var_sz = process_v(v)
      update(vcf_out, v, var_k, var_cat, 'SNP' if var_sz == 0 else 'INDEL', file_no, partitions, working_dict)
      if var_k is not None: ctr += 1

    if 0 <= ctr % 100000 <= 1:
      t1 = time.time()
      logger.debug('Read {} variants ({:0.2f} v/s). Working dict size is {}'.format(ctr, ctr/(t1 - t0), len(working_dict)))

  t1 = time.time()
  logger.debug(
    'Finished reading {} variants ({:0.2f} v/s). Residual working dict size is {}'.format(ctr, ctr / (t1 - t0), len(working_dict)))

  # At this point, the working_dict better only contain N -> FP and FP -> N
  for k, v in working_dict.items():
    v_to_use, cat_k = (v[1], 'N->FP') if v[0] is None else (v[0], 'FP->N')
    assert v_to_use['cat'] == 'FP', '{}: Not a FP. There is a bug in the code'.format(k)
    partitions[cat_k][v_to_use['type']] += 1
    write_out_v(vcf_out, v_to_use['v'], v_to_use['type'], cat_k)

  summary = summary_table_text(partitions)
  logger.debug(summary)

  summary_out.write(summary)


def summary_table_text(partitions):
  def line(k):
    return ['{}:\t\t{}\t{}'.format(k, partitions[k]['SNP'], partitions[k]['INDEL'])]

  res = []

  res += ['', 'Improved\tSNP\tINDEL', '-' * 20]
  for key in partition_list()[:4]:
    res += line(key)

  res += ['', 'Unchanged\tSNP\tINDEL', '-' * 20]
  for key in partition_list()[4:8]:
    res += line(key)

  res += ['', 'Regressed\tSNP\tINDEL', '-' * 20]
  for key in partition_list()[8:]:
    res += line(key)

  return '\n'.join(res)