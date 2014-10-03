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


class MyApp(wx.App):
  def OnInit(self):
    self.frame = MyFrame(None, title="Bitmaps")
    return True


  def in_lst(self, lst):
    self.frame.set_bmp(lst[0])
    print "in_lst"
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
    self.bitmap = wx.StaticBitmap(self.panel, bitmap = bmp_in)

if(__name__ == "__main__"):
  size_xy = 5
  data2d = flex.double(flex.grid(1, size_xy, size_xy),15)
  data2d[0, 1, 1] = 50

  for row in range(size_xy):
    for col in range(size_xy):
      data2d[0,row, col] += row * 2
      data2d[0,row, col] += col * 2

  mask2d = flex.int(flex.grid(1, size_xy, size_xy),size_xy)
  mask2d[0, 1, 1] = 5

  #testing wxbitmap_convert as a class
  a = wxbitmap_convert(data2d)
  print "calling obj", a.get_np()
  print a.__doc__


  datasize_xyd = flex.double(flex.grid(size_xy, size_xy, size_xy),15)
  datasize_xyd[1, 1, 1] = 50

  for frm in range(size_xy):
    for row in range(size_xy):
      for col in range(size_xy):
        datasize_xyd[frm, row, col] += (row * 2 + col * 2 + frm * 2)
        #datasize_xyd[0,row, col] += col * 2

  masksize_xyd = flex.int(flex.grid(size_xy, size_xy, size_xy),size_xy)
  masksize_xyd[0, 1, 1] = 5

  #testing wxbitmap_convert as a function
  print wxbitmap_convert(datasize_xyd).get_np()

  app = MyApp(redirect=False)
  app.in_lst(wxbitmap_convert(datasize_xyd).get_wxbitmap_lst())
  app.MainLoop()
