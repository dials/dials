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
import numpy as np
import wx.lib.scrolledpanel as scroll_pan

class show_3d(object):
  def __init__(self, data_xyz_in):
    app = show_3d_wx_app(redirect=False)
    app.in_lst(wxbitmap_convert(data_xyz_in).get_wxbitmap_lst())
    app.MainLoop()


class show_3d_wx_app(wx.App):
  def OnInit(self):
    self.frame = My_3d_flex_arr_frame(None, title="Bitmaps")
    return True


  def in_lst(self, lst):
    self.frame.ini_n_intro(lst)
    self.SetTopWindow(self.frame)
    self.frame.Show()


class multi_img_scrollable(scroll_pan.ScrolledPanel):
    def __init__(self, outer_frame, bmp_lst_in):
      super(multi_img_scrollable, self).__init__(outer_frame)

      self.set_scroll_content(bmp_lst_in)

      self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
      self.SetupScrolling()


    def set_scroll_content(self, bmp_lst_in):

      my_sizer = wx.BoxSizer(wx.HORIZONTAL)
      for bmp_lst in bmp_lst_in:
        local_bitmap = wx.StaticBitmap(self, bitmap = bmp_lst)
        my_sizer.Add(local_bitmap, 0, wx.LEFT | wx.ALL, 3)

      self.SetSizer(my_sizer)

      to_be_used_soon = '''
        for child in panel.GetChildren():
          child.Destroy()
      '''


    def OnMouseWheel(self, event):
        print "Weel event"
        print "event.GetWheelRotation() =", event.GetWheelRotation()


class buttons_panel(wx.Panel):
    def __init__(self, outer_frame):
        super(buttons_panel, self).__init__(outer_frame)
        self.p_frame  = outer_frame

        Hide_I_Button = wx.Button(self, label="Hide I")
        Hide_I_Button.Bind(wx.EVT_BUTTON, self.OnShwIBut)

        Show_I_Button = wx.Button(self, label="Show I")
        Show_I_Button.Bind(wx.EVT_BUTTON, self.OnHidIBut)

        my_sizer = wx.BoxSizer(wx.VERTICAL)
        my_sizer.Add(Show_I_Button, 0, wx.LEFT | wx.ALL,8)
        my_sizer.Add(Hide_I_Button, 0, wx.LEFT | wx.ALL,8)
        my_sizer.SetMinSize((90, 250))

        self.SetSizer(my_sizer)

    def OnShwIBut(self, event):
        print "OnShwIBut"
        self.p_frame.tst()

    def OnHidIBut(self, event):
        print "OnHidIBut"

class My_3d_flex_arr_frame(wx.Frame):
  def __init__(self, parent, id = wx.ID_ANY, title = "",
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style = wx.DEFAULT_FRAME_STYLE):
    super(My_3d_flex_arr_frame, self).__init__(parent, id, title,
                                  pos, size, style)

  def ini_n_intro(self, bmp_lst_in):

    # Attributes
    self.panel_01 = buttons_panel(self)
    self.panel_02 = multi_img_scrollable(self, bmp_lst_in)
    # Layout
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(self.panel_01, 0, wx.EXPAND)
    sizer.Add(self.panel_02, 1, wx.EXPAND)
    self.SetSizer(sizer)

    self.Show(True)
  def tst(self):
    print "testing"

if(__name__ == "__main__"):
  size_xy = 6

  data_xyz_flex = flex.double(flex.grid(size_xy, size_xy, size_xy),15)
  data_xyz_flex[1, 2, 2] = 15
  data_xyz_flex[2, 2, 2] = 20
  data_xyz_flex[3, 2, 2] = 25
  data_xyz_flex[4, 2, 2] = 20
  data_xyz_flex[5, 2, 2] = 15

  for frm in range(size_xy):
    for row in range(size_xy):
      for col in range(size_xy):
        data_xyz_flex[frm, row, col] += (row * 2 + col * 2 + frm * 2)


  show_3d(data_xyz_flex)
  print "Done"
