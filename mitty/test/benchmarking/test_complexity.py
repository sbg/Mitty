from mitty.benchmarking.complexity import ShannonEntropy, LinguisticComplexity, de_bruijn


def test_shannon():
  """complexity : Test shannon entropy for different k-mers"""
  for k in range(1, 6):
    se = ShannonEntropy(k=k)
    assert se.complexity(de_bruijn('ACTG', k)) == 2 * k


def test_lc():
  """complexity : Lingusitic complexity for de Bruijn sequences"""
  lc = LinguisticComplexity()
  for k in range(1, 5):
    assert lc.complexity(de_bruijn('ACTG', k)) == 1
