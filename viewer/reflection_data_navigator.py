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
from dials.viewer.viewer_utilities import build_np_img

class table_s_navigator(object):

  def __init__(self, table_in):
    self.table = table_in
    self.num_ref = len(self.table)
    if self.num_ref >= 1:
      self.row_pos = 0
    else:
      print "ERROR 0 reflections"
    self.z = 0
  def __call__(self):

    self.get_dat_bkg_msk()

    '''
    data2d = self.data_flex[self.z:self.z + 1, :, :]
    data2d.reshape(flex.grid(self.data_flex.all()[1:]))
    data2d_np = data2d.as_numpy_array()

    background2d = self.background_flex[self.z:self.z + 1, :, :]
    background2d.reshape(flex.grid(self.background_flex.all()[1:]))
    background2d_np = background2d.as_numpy_array()

    mask2d = self.mask_flex[self.z:self.z + 1, :, :]
    mask2d.reshape(flex.grid(self.mask_flex.all()[1:]))
    mask2d_np = mask2d.as_numpy_array()

    self.img_data = data2d_np
    self.img_background = background2d_np
    self.img_mask = mask2d_np
    return self.img_data, self.img_background, self.img_mask
    '''

    np_lst = []

    for post in range(3):
      z_from = self.z + post - 1
      z_to = self.z + post
      if z_from >= 0 and z_to <=self.depth:
        data2d = self.data_flex[z_from:z_to, :, :]
        data2d.reshape(flex.grid(self.data_flex.all()[1:]))
        data2d_np = data2d.as_numpy_array()
        np_lst.append(data2d_np)
      else:
        np_lst.append(build_np_img())

    '''
    self.img_data = np_lst[0]
    self.img_background = np_lst[1]
    self.img_mask = np_lst[2]

    return self.img_data, self.img_background, self.img_mask
    '''
    return np_lst


  def get_dat_bkg_msk(self):
    self.row = self.table[self.row_pos]
    self.data_flex = self.row['shoebox'].data
    self.background_flex = self.row['shoebox'].background
    self.mask_flex = self.row['shoebox'].mask
    self.box_lim = self.row['bbox']
    self.I_Max = flex.max(self.data_flex)
    self.depth = self.data_flex.all()[0]
    if self.depth <= 0:
      print "ERROR 0 depth"

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
  def Get_Max(self):
    return self.I_Max

  def Get_bbox(self):
    return self.box_lim

  '''
  def background(self):
    #print "from background(self)"
    return self.img_background

  def data(self):
    #print "from data(self)"
    return self.img_data

  def mask(self):
    #print "from mask(self)"
    return self.img_mask
  '''
