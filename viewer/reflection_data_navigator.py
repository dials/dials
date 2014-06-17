#
#  DIALS viewer
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

from dials.array_family import flex

class table_s_navigator(object):

  def __init__(self, table):

    self.table = table
    self.num_ref = len(table)
    if self.num_ref >= 1:
      row = table[52]
    else:
      print "ERROR 0 reflections"

    self.data_flex = row['shoebox'].data
    self.background_flex = row['shoebox'].background
    self.mask_flex = row['shoebox'].mask
    self.depth = self.data_flex.all()[0]
    print "depth of refl =", self.depth
    if self.depth >= 1:
      self.z = 0
    else:
      print "ERROR 0 depth"

  def next_slice(self):
    if self.z < self.depth:
      self.z += 1
      self.__call__()
    else:
      print "maximum depth reached"
  def Previous_slice(self):
    pass
  def next_Reflection(self):
    pass
  def Previous_Reflection(self):
    pass

  def __call__(self):

    mask2d = self.mask_flex[self.z:self.z + 1, :, :]
    mask2d.reshape(flex.grid(self.mask_flex.all()[1:]))
    mask2d_np = mask2d.as_numpy_array()

    data2d = self.data_flex[self.z:self.z + 1, :, :]
    data2d.reshape(flex.grid(self.data_flex.all()[1:]))
    data2d_np= data2d.as_numpy_array()

    background2d = self.background_flex[self.z:self.z + 1, :, :]
    background2d.reshape(flex.grid(self.background_flex.all()[1:]))
    background2d_np = background2d.as_numpy_array()

    self.img_background = background2d_np
    self.img_data = data2d_np
    self.img_mask = mask2d_np
    return self.img_background, self.img_data, self.img_mask

  def background(self):
    print "from background(self)"
    return self.img_background

  def data(self):
    print "from data(self)"
    return self.img_data

  def mask(self):
    print "from mask(self)"
    return self.img_mask