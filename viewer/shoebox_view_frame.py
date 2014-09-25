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

class ShoeboxView(wx.Frame):
  def __init__(self, parent, refl, orient = "portrait", id = wx.ID_ANY,
               title = "", pos = wx.DefaultPosition, size = wx.DefaultSize,
               style = wx.DEFAULT_FRAME_STYLE,
               name="ShoeboxView"):
    super(ShoeboxView, self).__init__(parent, id, title,
                                  pos, size, style, name)

    dat_flex = refl['shoebox'].data

    n_fr = dat_flex.all()[0]
    self.ImgLst = ImageListCtrl(self, orient)
    for fr in range(n_fr):
      data2d_flex = dat_flex[fr:fr + 1, :, :]
      data2d_flex.reshape(flex.grid(dat_flex.all()[1:]))
      data2d = data2d_flex.as_numpy_array()
      bitmap = GetBitmap_from_np_array(data2d)
      self.ImgLst.AppendBitmap(bitmap)


    self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.TstButton = wx.Button(self, label="Enbiggen")
    self.TstButton.Bind(wx.EVT_BUTTON, self.OnTstBut)
    self.top_sizer.Add(self.TstButton)

    self.Tst_01_Button = wx.Button(self, label="EnSmallen")
    self.Tst_01_Button.Bind(wx.EVT_BUTTON, self.OnTst_01_But)
    self.top_sizer.Add(self.Tst_01_Button)

    self.bot_sizer = wx.BoxSizer(wx.VERTICAL)
    self.bot_sizer.Add(self.top_sizer)
    self.bot_sizer.Add(self.ImgLst, 1, wx.EXPAND)

    self.SetSizer(self.bot_sizer)

  def OnTstBut(self, event):
    print "OnTstBut"

  def OnTst_01_But(self, event):
    print "OnTst_01_But"
