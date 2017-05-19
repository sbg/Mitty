from mitty.lib.cigars import cigarv2_v1


def test_cigar_conversion():
  """Cigars: Converting cigar V2 to V1"""
  test_cases = [
    ('33=1X79=1X26=1X109=', '250M'),
    ('1X26=1X123=1X82=1X15=', '250M'),
    ('89=10D161=', '89M10D161M'),
    ('99M1X', '100M'),
    ('10M10D1X9M', '10M10D10M')
  ]
  for tc in test_cases:
    assert cigarv2_v1(tc[0]) == tc[1], cigarv2_v1(tc[0])