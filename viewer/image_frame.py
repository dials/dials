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

from dials.viewer.img_utilities import GetBitmap_from_np_array, ImageListCtrl

from dials.array_family import flex

class ImageFrame(wx.Frame):
  def __init__(self, parent, refl, id=wx.ID_ANY, title="",
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style=wx.DEFAULT_FRAME_STYLE,
               name="ImageFrame"):
    super(ImageFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)

    dat_flex = refl['shoebox'].data


    n_fr = dat_flex.all()[0]
    self.ImgLst = ImageListCtrl(self)
    for fr in range(n_fr):
      data2d_flex = dat_flex[fr:fr + 1, :, :]
      data2d_flex.reshape(flex.grid(dat_flex.all()[1:]))
      data2d = data2d_flex.as_numpy_array()
      bitmap = GetBitmap_from_np_array(data2d)
      self.ImgLst.AppendBitmap(bitmap)

    #self.panel = wx.Panel(self)
    #self.bitmap = wx.StaticBitmap(self.panel, bitmap = bitmap)
    self.mid_sizer = wx.BoxSizer(wx.VERTICAL)
    self.mid_sizer.Add(self.ImgLst, 1, wx.EXPAND)
    self.SetSizer(self.mid_sizer)


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
