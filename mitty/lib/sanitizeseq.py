"""Remove all non ACTG bases from a sequence"""
sanitable = str.maketrans('ATCGNRYSWKMBDHV', 'ATCGNNNNNNNNNNN')


def sanitize(seq):
  return seq.translate(sanitable)