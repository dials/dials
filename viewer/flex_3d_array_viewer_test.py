#!/usr/bin/env python
#
#  flex_3d_array_viewer_test.py
#
#  test for slice_viewer.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
from dials.array_family import flex
from dials.viewer.slice_viewer import show_3d
from dials.algorithms.shoebox import MaskCode

attempt_to_fix_a_crash = '''
try:
  import scipy.linalg # import dependency
except ImportError, e:
  pass
'''

if(__name__ == "__main__"):
  lst_flex = []
  lst_flex_norm = []

  for size_xyz in range(12, 7, -1):

    data_xyz_flex = flex.double(flex.grid(size_xyz, size_xyz, size_xyz),15)
    data_flex_norm = flex.double(flex.grid(size_xyz, size_xyz, size_xyz),15)
    data_flex_mask = flex.int(flex.grid(size_xyz, size_xyz, size_xyz),0)

    #data_xyz_flex[1, 2, 2] = 35 + size_xyz * 5
    #data_xyz_flex[2, 2, 2] = 40 + size_xyz * 5

    tot = 0.0
    for frm in range(size_xyz):
      for row in range(size_xyz):
        for col in range(size_xyz):
          data_xyz_flex[frm, row, col] += (row * 2 + col * 2 + frm * 2)
          tot += data_xyz_flex[frm, row, col]
          if( row > 1 and row < size_xyz - 2 and col > 1 and col < size_xyz - 2 ):
            data_flex_mask[frm, row, col] = MaskCode.Foreground


            different_mask_values = '''
            MaskCode.Valid           =  "\\\\\\"
            MaskCode.Foreground      =  "//////"
            MaskCode.Background      =  "||||||"
            MaskCode.BackgroundUsed  =  "------"
            '''


    for frm in range(size_xyz):
      for row in range(size_xyz):
        for col in range(size_xyz):
          data_flex_norm[frm, row, col] += data_xyz_flex[frm, row, col] / tot


    lst_flex.append(data_xyz_flex)
    lst_flex_norm.append(data_flex_norm)


  show_3d(data_xyz_flex, data_flex_mask)
  show_3d(lst_flex)
  show_3d(lst_flex_norm)
