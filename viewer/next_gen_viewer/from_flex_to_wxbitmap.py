#!/usr/bin/env python
#
# dials.reflection_viewer.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.array_family import flex

from bitmap_from_numpy_w_matplotlib import GetBitmap_from_np_array
class wxbitmap_convert(object):
  '''
  The main duty of this class is to convert from
  a 3D flex array to a list of WxBitmaps
  '''
  def __init__(self, data_in):
    self.lst_np = []

    try:
      self.depth = data_in.all()[0]
    except:
      print "error in entered data"
      self.depth = 0

    for z_from in range(self.depth):
      data2d = data_in[z_from:z_from + 1, :, :]
      data2d.reshape(flex.grid(data_in.all()[1:]))
      data2d_np = data2d.as_numpy_array()
      self.lst_np.append(data2d_np)


  def get_np_lst(self):
    return self.lst_np


  def get_wxbitmap_lst(self):
    lst_img = []
    for np_elm in self.lst_np:
      wxbmp = GetBitmap_from_np_array(np_elm)
      lst_img.append(wxbmp)
    return lst_img
