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

from bitmap_from_numpy_w_matplotlib import wxbmp_from_np_array
class wxbitmap_convert(object):
  '''
  The main duty of this class is to convert from
  a 3D flex array to a list of WxBitmaps
  '''
  def __init__(self, data_in):

    try:
      self.depth = data_in.all()[0]
    except:
      print "error in entered data"
      self.depth = 0

    data2d = data_in[0:1, :, :]
    data2d_np = data2d.as_numpy_array()
    self.np_3d_block = data2d_np


  def get_np(self):
    return self.np_3d_block


  def get_wxbitmap_lst(self):
    lst_img = []
    local_bmp = wxbmp_from_np_array()
    wxbmp = local_bmp.get_bmp(self.np_3d_block)
    lst_img.append(wxbmp)
    return lst_img
