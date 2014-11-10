#
#  DIALS reflection table navigator
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."
from __future__ import division
from dials.array_family import flex

class table_s_navigator(object):

  def __init__(self, table_in):
    self.table = table_in
    self.num_ref = len(self.table)
    if self.num_ref >= 1:
      self.row_pos = 0
    else:
      #print "ERROR 0 reflections"
      pass
    self.z = 0
    print "Num of Refls = ", self.num_ref

  def __call__(self, opt = 0):
    self.get_dat_bkg_msk()
    np_lst = []
    title_lst = []
    if( opt == 0 ):
      try:
        data2d = self.data_flex[self.z:self.z + 1, :, :]
        data2d.reshape(flex.grid(self.data_flex.all()[1:]))
        np_lst.append(data2d.as_numpy_array())
        title_lst.append("shoebox_data[" + str(self.z) + ":" + str(self.z+1) + "]")
      except:
        np_lst.append(None)
        title_lst.append("No data")
      try:
        background2d = self.background_flex[self.z:self.z + 1, :, :]
        background2d.reshape(flex.grid(self.background_flex.all()[1:]))
        np_lst.append(background2d.as_numpy_array())
        title_lst.append("shoebox_background[" + str(self.z) + ":" + str(self.z+1) + "]")
      except:
        np_lst.append(None)
        title_lst.append("No background")

      try:
        mask2d = self.mask_flex[self.z:self.z + 1, :, :]
        mask2d.reshape(flex.grid(self.mask_flex.all()[1:]))
        np_lst.append(mask2d.as_numpy_array())
        title_lst.append("shoebox_mask[" + str(self.z) + ":" + str(self.z+1) + "]")
      except:
        np_lst.append(None)
        title_lst.append("No mask")

    elif( opt == 1 or opt == 2 or opt == 3 ):

      for post in range(3):
        z_from = self.z + post - 1
        z_to = self.z + post
        if z_from >= 0 and z_to <=self.depth:
          if( opt == 1 ):
            try:
              data2d = self.data_flex[z_from:z_to, :, :]
              title_lst.append("shoebox_data[" + str(z_from) + ":" + str(z_to) + "]")
            except:
              np_lst.append(None)
          elif( opt == 2 ):
            try:
              data2d = self.background_flex[z_from:z_to, :, :]
              title_lst.append("shoebox_background[" + str(z_from) + ":" + str(z_to) + "]")
            except:
              np_lst.append(None)
          elif( opt == 3 ):
            try:
              data2d = self.mask_flex[z_from:z_to, :, :]
              title_lst.append("shoebox_mask[" + str(z_from) + ":" + str(z_to) + "]")
            except:
              np_lst.append(None)
          try:
            data2d.reshape(flex.grid(self.data_flex.all()[1:]))
            data2d_np = data2d.as_numpy_array()
            np_lst.append(data2d_np)
          except:
            pass


        else:
          np_lst.append(None)
          title_lst.append(None)
    else:
      #print "wrong option"
      np_lst = [None, None, None]
      title_lst = [None, None, None]
    return np_lst, title_lst

  def get_dat_bkg_msk(self):
    try:
      table_row = self.table[self.row_pos]
    except:
      #print "No Table to read from"
      pass

    try:
      self.data_flex = table_row['shoebox'].data
    except:
      self.data_flex = None
      #print "No shoebox IMG data"

    try:
      self.background_flex = table_row['shoebox'].background
    except:
      self.background_flex = None
      #print "No background data"
      #self.background_flex = None


    try:
      self.hkl_data = table_row['miller_index']
    except:
      self.hkl_data = None
      #print "No HKL data"

    try:
      self.mask_flex = table_row['shoebox'].mask
    except:
      self.mask_flex = None
      #print "No mask data"

    try:
      self.r_bbox = table_row['bbox']
    except:
      self.r_bbox = None
      #print "No bbox data"

    try:
      self.xyzcal_px = list(table_row['xyzcal.px'])
      #'''
      self.xyzcal_px[0] = self.xyzcal_px[0] - self.r_bbox[0]
      self.xyzcal_px[1] = self.xyzcal_px[1] - self.r_bbox[2]
      #'''
    except:
      self.xyzcal_px = None
      #print "No xyzcal.px data"

    try:
      self.depth = self.data_flex.all()[0]
    except:
      self.depth = 0

    if self.depth <= 0:
      #print "ERROR 0 depth"
      pass
    return table_row

  def next_slice(self):
    if self.z < self.depth - 1:
      self.z += 1
    else:
      #print "maximum depth reached"
      pass

  def Previous_slice(self):
    if self.z > 0:
      self.z -= 1
    else:
      #print "depth 0 reached"
      pass

  def next_Reflection(self):
    if self.row_pos < self.num_ref - 1:
      self.row_pos += 1
      self.z = 0
    else:
      #print "last reflection reached"
      pass

  def Previous_Reflection(self):
    if self.row_pos > 0:
      self.row_pos -= 1
      self.z = 0
    else:
      #print "first reflection reached"
      pass

  def Jump_to_Reflection(self, Ref_Num):
    if Ref_Num > self.num_ref - 1 or Ref_Num < 0:
      #print "There is no reflection #", Ref_Num
      pass
    else:
      self.row_pos = Ref_Num
      self.z = 0

  def Get_Max(self, opt = 0):
    if( self.data_flex !=None ):
      if(opt != 3):
        self.I_Max = flex.max(self.data_flex)
      else:
        self.I_Max = flex.max(self.mask_flex) * 2
    else:
      self.I_Max = -1
    return self.I_Max

  def Get_xyz(self):
    return self.xyzcal_px
  def Get_hkl(self):
    return self.hkl_data
  def Get_ref_num(self):
    return self.row_pos

  def Get_bbox(self):
    return self.r_bbox

