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
from from_flex_to_wxbitmap import wxbitmap_convert
import wx

class show_3d(object):
  def __init__(self, data_xyz_in):
    app = show_3d_wx_app(redirect=False)
    app.in_lst(wxbitmap_convert(data_xyz_in).get_wxbitmap_lst())
    app.MainLoop()


class show_3d_wx_app(wx.App):
  def OnInit(self):
    self.frame = MyFrame(None, title="Bitmaps")
    return True


  def in_lst(self, lst):
    self.frame.set_bmp(lst)
    self.SetTopWindow(self.frame)
    self.frame.Show()


class MyFrame(wx.Frame):
  def __init__(self, parent, id = wx.ID_ANY, title = "",
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style = wx.DEFAULT_FRAME_STYLE,
               name = "MyFrame"):
    super(MyFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)
    # Attributes
    self.panel = wx.Panel(self)

  def set_bmp(self, bmp_in):
    self.img_sizer = wx.BoxSizer(wx.HORIZONTAL)

    for bmp_lst in bmp_in:
      local_bitmap = wx.StaticBitmap(self.panel, bitmap = bmp_lst)
      self.img_sizer.Add(local_bitmap, 0, wx.LEFT | wx.ALL,3)


    self.panel.SetSizerAndFit(self.img_sizer)

if(__name__ == "__main__"):
  size_xy = 5

  data_xyz_flex = flex.double(flex.grid(size_xy, size_xy, size_xy),15)
  data_xyz_flex[2, 1, 1] = 50

  for frm in range(size_xy):
    for row in range(size_xy):
      for col in range(size_xy):
        data_xyz_flex[frm, row, col] += (row * 2 + col * 2 + frm * 2)


  show_3d(data_xyz_flex)
  print "Done"
