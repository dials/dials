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
  lst_flex_norm = []
  for size_xyz in range(3,8):

    data_xyz_flex = flex.double(flex.grid(size_xyz, size_xyz, size_xyz),15)
    data_flex_norm = flex.double(flex.grid(size_xyz, size_xyz, size_xyz),15)
    #data_xyz_flex[1, 2, 2] = 35 + size_xyz * 5
    #data_xyz_flex[2, 2, 2] = 40 + size_xyz * 5

    tot = 0.0
    for frm in range(size_xyz):
      for row in range(size_xyz):
        for col in range(size_xyz):
          data_xyz_flex[frm, row, col] += (row * 2 + col * 2 + frm * 2)
          tot += data_xyz_flex[frm, row, col]

    for frm in range(size_xyz):
      for row in range(size_xyz):
        for col in range(size_xyz):
          data_flex_norm[frm, row, col] += data_xyz_flex[frm, row, col] / tot


    lst_flex.append(data_xyz_flex)
    lst_flex_norm.append(data_flex_norm)


  show_3d(data_xyz_flex)

  show_3d(lst_flex)
  show_3d(lst_flex_norm)
