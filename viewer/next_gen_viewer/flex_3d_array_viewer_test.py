#!/usr/bin/env python
#
#  flex_3d_array_viewer_test.py
#
#  test for multi_3D_slice_viewer.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from dials.array_family import flex
from dials.viewer.next_gen_viewer.multi_3D_slice_viewer import show_3d
if(__name__ == "__main__"):
  lst_flex = []
  for size_xy in range(3,7):

    data_xyz_flex = flex.double(flex.grid(size_xy, size_xy, size_xy),15)
    data_xyz_flex[1, 2, 2] = 35
    data_xyz_flex[2, 2, 2] = 40

    for frm in range(size_xy):
      for row in range(size_xy):
        for col in range(size_xy):
          data_xyz_flex[frm, row, col] += (row * 2 + col * 2 + frm * 2)

    lst_flex.append(data_xyz_flex)

  show_3d(data_xyz_flex)

  show_3d(lst_flex)

