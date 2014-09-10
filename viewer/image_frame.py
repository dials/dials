#
#  DIALS Image Utilities
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luiso and James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."
#
from __future__ import division
import wx
#from dials.scratch.luiso_s.wx_toys.bitmap_from_numpy_w_matplotlib_well_done \
#     import GetBitmap_from_np_array

from dials.viewer.img_utilities import GetBitmap_from_np_array

from dials.array_family import flex

class ImageFrame(wx.Frame):
  def __init__(self, parent, refl, id=wx.ID_ANY, title="",
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style=wx.DEFAULT_FRAME_STYLE,
               name="ImageFrame"):
    super(ImageFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)
    print refl

    dat_flex = refl['shoebox'].data

    data2d_flex = dat_flex[0:1, :, :]
    data2d_flex.reshape(flex.grid(dat_flex.all()[1:]))

    # Attributes
    data2d = data2d_flex.as_numpy_array()

    print "data2d ="
    print data2d

    bitmap = GetBitmap_from_np_array(data2d)

    self.panel = wx.Panel(self)
    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)

  '''
  @property
  def image(self):
    return self._image

  @image.setter
  def image(self, value):

      #Value is a 2d flex array


    #return self._image = value
    pass
  '''
