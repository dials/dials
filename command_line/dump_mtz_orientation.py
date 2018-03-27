#!/usr/bin/env python
# LIBTBX_SET_DISPATCHER_NAME dev.dials.dump_mtz_orientation
#
# dials.command_line.dump_mtz_orientation
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
from __future__ import print_function

def dump_mtz_orientation(mtz_file):
  from iotbx import mtz
  from scitbx import matrix
  from scitbx.math.euler_angles import xyz_angles

  m = mtz.object(mtz_file)
  for b in m.batches():
    rxyz = tuple(xyz_angles(matrix.sqr(b.umat())))
    print(b.num(), '%7.4f %7.4f %7.4f' % rxyz)

if __name__ == '__main__':
  import sys
  for mtz_file in sys.argv[1:]:
    dump_mtz_orientation(mtz_file)
