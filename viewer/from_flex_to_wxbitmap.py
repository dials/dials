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
import numpy as np

from bitmap_from_numpy_w_matplotlib import wxbmp_from_np_array
class wxbitmap_convert(object):
  '''
  The main duty of this class is to convert
  from
  a 3D flex array or a list of flex arrays
  to
  a list of WxBitmaps
  '''
  def __init__(self, data_in_n1, data_in_n2 = None):

    if( data_in_n2 == None ):
      self.lst_3d_mask = None
      print "No double list given"

      if type(data_in_n1) is list:
        print "is a single list"
        self.lst_3d_data = []
        for lst_memb in data_in_n1:
          self.lst_3d_data.append(lst_memb.as_numpy_array())

      else:
        print "Got flex array"
        self.lst_3d_data = []
        self.lst_3d_data.append(data_in_n1.as_numpy_array())


    else:

      print "Got two arguments"
      if( type(data_in_n1) is list and type(data_in_n2) is list):
        print "Got two lists"
        if( len(data_in_n1) == len(data_in_n2) ):
          self.lst_3d_data = []
          self.lst_3d_mask = []
          for lst_pos in range(len(data_in_n1)):
            lst_memb1 = data_in_n1[lst_pos]
            self.lst_3d_data.append(lst_memb1.as_numpy_array())
            lst_memb2 = data_in_n2[lst_pos]
            self.lst_3d_mask.append(lst_memb2.as_numpy_array())

        else:
          print "the two lists do NOT have the same size"

      elif( type(data_in_n1) is not list and type(data_in_n2) is not list ):
        print "Got two blocks"
        self.lst_3d_data = []
        self.lst_3d_mask = []
        self.lst_3d_data.append(data_in_n1.as_numpy_array())

        self.lst_3d_mask.append(  data_in_n2.as_numpy_array() )

        old_testing_stuff = '''
        print "nump_arr =", nump_arr
        zmax = nump_arr.shape[0]
        xmax = nump_arr.shape[1]
        ymax = nump_arr.shape[2]
        print "xmax, ymax, zmax =", xmax, ymax, zmax
        print "len(nump_arr) =", len(nump_arr)
        nump_3d = np.zeros( (zmax, xmax, ymax), 'int')

        nump_3d[:,:,:] = nump_arr[:,:,:]
        self.lst_3d_mask.append(  nump_3d )
        '''

      else:
        print "Got mixture of different type of data"


  def get_np(self):
    return self.lst_3d_data


  def get_wxbitmap_lst(self, show_nums = True, scale = 1.0):
    self.local_bmp = wxbmp_from_np_array(self.lst_3d_data, show_nums, self.lst_3d_mask)
    return self.scaling(scale)


  def scaling(self, scale = 1.0):
    lst_img = self.local_bmp.bmp_lst_scaled(scale)


    return lst_img
