"""Functions to help parse an EVCF as we want"""
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
  sz = len(alt_dict[cat](t, q)) - len(v.ref) if cat is not None and cat is not 'skip' else None
  return cat, sz