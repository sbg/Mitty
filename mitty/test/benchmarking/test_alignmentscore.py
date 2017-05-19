from mitty.benchmarking.alignmentscore import find_first_common_reference_matching_base_positions, pysam


def test_ref_matching_base():
  """d-err: Find first reference matching base positions"""
  test_cases = [
    ((100, '100M'), (130, '30S70M')),
    ((100, '10I90M'), (100, '10S90M')),
    ((100, '10I10M10I70M'), (110, '30S70M')),
    ((100, '10M1000D90M'), (1110, '10S90M')),
    ((100, '10M1000D90M'), (1120, '20S90M')),
    ((100, '10M1000D80M10I'), (1120, '20S80M10S')),
    ((1000686, '96S154M'), (1000590, '250M')),
  ]
  for tc1, tc2 in test_cases:
    r1, r2 = pysam.AlignedSegment(), pysam.AlignedSegment()
    r1.pos, r1.cigarstring = tc1[0], tc1[1]
    r2.pos, r2.cigarstring = tc2[0], tc2[1]

    p1, p2 = find_first_common_reference_matching_base_positions(r1, r2)

    assert p1 == p2 and p1 is not None, (tc1, tc2, p1, p2)