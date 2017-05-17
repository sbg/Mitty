import re

pattern = re.compile('([MIDNSHPX=])')


def cigarv2_v1(cigar_v2):
  """Convert a V2 cigar_v1 string to V1

  :param cigarstring:
  :return:
  """
  values = pattern.split(cigar_v2.replace('=', 'M').replace('X', 'M'))
  cigar_v1 = []
  last_op, op_cnt = values[1], 0
  for op in zip(values[::2], values[1::2]):
    if op[1] == last_op:
      op_cnt += int(op[0])
    else:
      cigar_v1 += [str(op_cnt), last_op]
      last_op, op_cnt = op[1], int(op[0])
  cigar_v1 += [str(op_cnt), last_op]
  return ''.join(cigar_v1)