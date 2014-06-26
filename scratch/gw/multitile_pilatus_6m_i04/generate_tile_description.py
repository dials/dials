from __future__ import division
def import_offsets():

  offsets = { }

  for record in open('offsets.txt'):
    tokens = record.split()
    if not tokens:
      continue
    row, column = int(tokens[0]), int(tokens[1])
    x, y, r = map(float, tokens[2:])
    offsets[(row, column)] = (x, y, r)

  return offsets

def print_offsets_in_cif(offsets):

  # first print the axis definitions

  for row, column in sorted(offsets):
    element_x = 'MOD_R%02d_C%02d_X' % (row, column)
    element_y = 'MOD_R%02d_C%02d_Y' % (row, column)
    element_r = 'MOD_R%02d_C%02d_R' % (row, column)

    offset_x = 0.172 * column * (487 + 7)
    offset_y = -0.172 * row * (195 + 17)

    print '%s translation detector DETECTOR_X    1 0 0 %8.2f %8.2f   0.00' % \
          (element_x, offset_x, offset_y)
    print '%s translation detector %s 0 1 0     0.00     0.00   0.00' % \
          (element_y, element_x)
    print '%s rotation    detector %s 0 0 1     0.00     0.00   0.00' % \
          (element_r, element_y)

  # now apply the shifts

  for row, column in sorted(offsets):

    offset = offsets[(row, column)]

    element_x = 'MOD_R%02d_C%02d_X' % (row, column)
    element_y = 'MOD_R%02d_C%02d_Y' % (row, column)
    element_r = 'MOD_R%02d_C%02d_R' % (row, column)

    import math
    r2d = 180.0 / math.pi

    print 'FRAME1 %s  0.00    %9.6f' % (element_x, offset[0] * 0.172)
    print 'FRAME1 %s  0.00    %9.6f' % (element_y, offset[1] * 0.172)
    print 'FRAME1 %s %9.6f 0.00' % (element_r, offset[2] * r2d)

def print_offsets_in_cif2(offsets):

  # first print the axis definitions

  for row, column in sorted(offsets):
    element = 'MOD_R%02d_C%02d' % (row, column)
    print ' %s %%(detector_id)s' % element

  # this uses (1) fortran ordering i.e. fast, slow not C slow, fast and also
  # is inclusive i.e. counts from 1 not 0?

  for row, column in sorted(offsets):
    offset_x = column * (487 + 7)
    offset_y = row * (195 + 17)
    element = 'MOD_R%02d_C%02d' % (row, column)
    array = 'ARRAY1(%d:%d,%d:%d)' % (offset_x + 1, offset_x + 487,
                                     offset_y + 1, offset_y + 195)
    print ' FRAME1 %s %-27s 1' % \
          (element, array)

  for row, column in sorted(offsets):
    offset_x = column * (487 + 7)
    offset_y = row * (195 + 17)
    array = 'ARRAY1(%d:%d,%d:%d)' % (offset_x + 1, offset_x + 487,
                                     offset_y + 1, offset_y + 195)
    element_x = 'MOD_R%02d_C%02d_X' % (row, column)
    element_y = 'MOD_R%02d_C%02d_Y' % (row, column)
    print ' %-27s 1 487 1 increasing %s' % (array, element_x)
    print ' %-27s 2 195 2 increasing %s' % (array, element_y)

  for row, column in sorted(offsets):
    element_x = 'MOD_R%02d_C%02d_X' % (row, column)
    element_y = 'MOD_R%02d_C%02d_Y' % (row, column)
    print ' %s %s 0.0  0.172' % (element_x, element_x)
    print ' %s %s 0.0 -0.172' % (element_y, element_y)

  for row, column in sorted(offsets):
    offset_x = column * (487 + 7)
    offset_y = row * (195 + 17)
    array = 'ARRAY1(%d:%d,%d:%d)' % (offset_x + 1, offset_x + 487,
                                     offset_y + 1, offset_y + 195)
    print ' %-27s 1 0.000172' % array
    print ' %-27s 2 0.000172' % array

if __name__ == '__main__':
  offsets = import_offsets()
  print_offsets_in_cif2(offsets)
