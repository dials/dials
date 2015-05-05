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
#from dials_viewer_ext import gen_img, rgb_img
from dials_viewer_ext import rgb_img, gen_str_tst
if(__name__ == "__main__"):
  tmp_not_needed = '''
  size_xyz = 7

  arr_2d = flex.double(flex.grid(size_xyz, size_xyz), 00)

  tot = 0.0
  for row in range(size_xyz):
    for col in range(size_xyz):
      arr_2d[row, col] += (row * 2 + col * 2)
      tot += arr_2d[row, col]


  print "arr_2d.as_numpy_array() = \n", arr_2d.as_numpy_array()
  #print "gen_img(arr_2d).as_numpy_array() = \n", gen_img(arr_2d).as_numpy_array()


  wx_bmp_arr = rgb_img()
  print wx_bmp_arr.gen_bmp(arr_2d).as_numpy_array()
  '''

  print "testing new Cpp code"
  size_lng = 25
  arr_1d = flex.double(flex.grid(size_lng, 1), 00)
  for i in xrange(size_lng):
    arr_1d[i, 0] = 99.0 * float(i * i * i) + 0.986764543
  gen_str_tst(arr_1d)
  #print "gen_str_tst(arr_1d).as_numpy_array() = \n", gen_str_tst(arr_1d).as_numpy_array()
