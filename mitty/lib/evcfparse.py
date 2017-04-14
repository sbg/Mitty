"""Functions to help parse an EVCF as we want


Helpful docs:

https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md


"""
cat_dict = {
  'TPTP': 'TP',   # These are all seen in combination with alternate representations
  'TP.': 'TP',
  'TPFP': 'TP',
  '.FP': 'FP',
  'FN.': 'FN',
  'FNFP': 'GT',
  '.TP': 'skip',  # These are all combinations that are scored TP/nocall later
  'FNTP': 'skip', # because they are different representations
  'N.': 'skip',   # Unclear about
  '.N': 'skip',
  'NN': 'skip'
}

alt_dict = {
  'TP': lambda t, q: get_alt(t),
  'FN': lambda t, q: get_alt(t),
  'GT': lambda t, q: get_alt(t),
  'FP': lambda t, q: get_alt(q),
}


def get_alt(s):
  for n, alt in zip(s.allele_indices, s.alleles):
    if n:
      return alt


def parse_line(v):
  """Given a variant (line) from the VCF categorize it as a TP, FP, FN or GT

  :param v:
  :return:   (category,
              size)
  """
  t, q = v.samples['TRUTH'], v.samples['QUERY']
  cat = cat_dict.get(t['BD'] + q['BD'], None)
  if cat is not None and cat is not 'skip':
    alt = alt_dict[cat](t, q)
    sz = len(alt) - len(v.ref)
  else:
    alt, sz = None, None

  return cat, sz, alt