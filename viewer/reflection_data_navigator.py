#
#  DIALS reflection table navigator
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

from dials.array_family import flex


class table_s_navigator(object):

  def __init__(self, table_in):
    self.table = table_in
    self.num_ref = len(self.table)
    if self.num_ref >= 1:
      self.row_pos = 0
    else:
      print "ERROR 0 reflections"
    self.z = 0

  def __call__(self, opt = 0):
    self.get_dat_bkg_msk()
    np_lst = []
    if( opt == 0 ):
      data2d = self.data_flex[self.z:self.z + 1, :, :]
      data2d.reshape(flex.grid(self.data_flex.all()[1:]))
      np_lst.append(data2d.as_numpy_array())

      background2d = self.background_flex[self.z:self.z + 1, :, :]
      background2d.reshape(flex.grid(self.background_flex.all()[1:]))
      np_lst.append(background2d.as_numpy_array())

      mask2d = self.mask_flex[self.z:self.z + 1, :, :]
      mask2d.reshape(flex.grid(self.mask_flex.all()[1:]))
      np_lst.append(mask2d.as_numpy_array())

    elif( opt == 1 or opt == 2 or opt == 3 ):

      for post in range(3):
        z_from = self.z + post - 1
        z_to = self.z + post
        if z_from >= 0 and z_to <=self.depth:
          if( opt == 1 ):
            data2d = self.data_flex[z_from:z_to, :, :]
          elif( opt == 2 ):
            data2d = self.background_flex[z_from:z_to, :, :]
          elif( opt == 3 ):
            data2d = self.mask_flex[z_from:z_to, :, :]
          data2d.reshape(flex.grid(self.data_flex.all()[1:]))
          data2d_np = data2d.as_numpy_array()
          np_lst.append(data2d_np)
        else:
          np_lst.append(None)
    else:
      print "wrong option"
      np_lst = [None, None, None]
    return np_lst

  def get_dat_bkg_msk(self):
    try:
      table_row = self.table[self.row_pos]
    except:
      print "No Table to read from"

    try:
      self.data_flex = table_row['shoebox'].data
    except:
      print "No shoebox IMG data"

    try:
      self.background_flex = table_row['shoebox'].background
    except:
      print "No background data"
      #self.background_flex = None

    try:
      self.mask_flex = table_row['shoebox'].mask
    except:
      print "No mask data"

    try:
      self.box_lim = table_row['bbox']
    except:
      print "No bbox data"

    try:
      self.calc_pos = table_row['xyzcal.px']
    except:
      print "No xyzcal.px data"
      self.calc_pos = None

    self.depth = self.data_flex.all()[0]
    if self.depth <= 0:
      print "ERROR 0 depth"
    return table_row

  def next_slice(self):
    if self.z < self.depth - 1:
      self.z += 1
    else:
      print "maximum depth reached"

  def Previous_slice(self):
    if self.z > 0:
      self.z -= 1
    else:
      print "depth 0 reached"

  def next_Reflection(self):
    if self.row_pos < self.num_ref - 1:
      self.row_pos += 1
      self.z = 0
    else:
      print "last reflection reached"

  def Previous_Reflection(self):
    if self.row_pos > 0:
      self.row_pos -= 1
      self.z = 0
    else:
      print "first reflection reached"

  def Jump_to_Reflection(self, Ref_Num):
    if Ref_Num > self.num_ref - 1 or Ref_Num < 0:
      print "There is no reflection #", Ref_Num
    else:
      self.row_pos = Ref_Num
      self.z = 0

  def Get_Max(self, opt = 0):
    if(opt != 3):
      self.I_Max = flex.max(self.data_flex)
    else:
      self.I_Max = flex.max(self.mask_flex) * 2
    return self.I_Max

  def Get_xyz(self):
    return self.calc_pos

  def Get_bbox(self):
    return self.box_lim

