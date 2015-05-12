#!/usr/bin/env python
#
#  boost_tst.py
#
#  test of generation of RGB images in C++
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.array_family import flex

import boost.python
from dials_viewer_ext import gen_font_img

if(__name__ == "__main__"):
  lst_flex = []
  lst_flex_norm = []

  depth = 16

  arr_1d = flex.double(flex.grid(depth), 00)

  for i in range(depth):
    arr_1d[i] = float(i)


  print "arr_1d.as_numpy_array() = \n", arr_1d.as_numpy_array()

  block_3d = gen_font_img(arr_1d).as_numpy_array()

  for i in range(depth):
    tmp_2d_arr = block_3d[:,:,i]
    print "arr 2d (", i, ") =\n", tmp_2d_arr

