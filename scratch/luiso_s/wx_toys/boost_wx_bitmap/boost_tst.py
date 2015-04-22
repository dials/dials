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
from dials_viewer_ext import gen_img

if(__name__ == "__main__"):
  lst_flex = []
  lst_flex_norm = []

  size_xyz = 7

  arr_2d = flex.double(flex.grid(size_xyz, size_xyz), 00)

  tot = 0.0
  for row in range(size_xyz):
    for col in range(size_xyz):
      arr_2d[row, col] += (row * 2 + col * 2)
      tot += arr_2d[row, col]


print "arr_2d.as_numpy_array() = \n", arr_2d.as_numpy_array()
print "gen_img(arr_2d).as_numpy_array() = \n", gen_img(arr_2d).as_numpy_array()

