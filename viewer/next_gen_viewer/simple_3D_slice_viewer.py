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
    app.in_lst(data_xyz_in)
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
      self.p_frame  = outer_frame
      self.local_bmp_lst = bmp_lst_in

      self.set_scroll_content(self.local_bmp_lst)

      self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
      self.SetupScrolling()


    def set_scroll_content(self, bmp_lst_in):

      self.my_img_sizer = wx.BoxSizer(wx.HORIZONTAL)
      for bmp_lst in bmp_lst_in:
        local_bitmap = wx.StaticBitmap(self, bitmap = bmp_lst)
        self.my_img_sizer.Add(local_bitmap, 0, wx.LEFT | wx.ALL, 3)

      self.SetSizer(self.my_img_sizer)


    def OnMouseWheel(self, event):
      rot_sn = event.GetWheelRotation()
      self.p_frame._to_re_zoom(rot_sn)


    def img_refresh(self, bmp_lst_new):
      self.local_bmp_lst = bmp_lst_new
      for child in self.GetChildren():
        child.Destroy()
      self.set_scroll_content(self.local_bmp_lst)
      self.Layout()
      self.p_frame.Layout()
      self.Refresh()


class buttons_panel(wx.Panel):
    def __init__(self, outer_frame):
        super(buttons_panel, self).__init__(outer_frame)
        self.p_frame  = outer_frame

        Hide_I_Button = wx.Button(self, label="Hide I")
        Hide_I_Button.Bind(wx.EVT_BUTTON, self.OnHidIBut)
        Show_I_Button = wx.Button(self, label="Show I")
        Show_I_Button.Bind(wx.EVT_BUTTON, self.OnShwIBut)

        self.my_sizer = wx.BoxSizer(wx.VERTICAL)
        self.my_sizer.Add(Show_I_Button, 0, wx.LEFT | wx.ALL,8)
        self.my_sizer.Add(Hide_I_Button, 0, wx.LEFT | wx.ALL,8)
        self.my_sizer.SetMinSize((90, 250))

        self.SetSizer(self.my_sizer)

    def OnHidIBut(self, event):
        self.p_frame._to_hide_nums()

    def OnShwIBut(self, event):
        self.p_frame._to_show_nums()

class My_3d_flex_arr_frame(wx.Frame):
  def __init__(self, parent, id = wx.ID_ANY, title = "",
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style = wx.DEFAULT_FRAME_STYLE):
    super(My_3d_flex_arr_frame, self).__init__(parent, id, title,
                                  pos, size, style)
    self.show_nums = True

  def ini_n_intro(self, flex_arr_in):
    self.flex_arr = flex_arr_in
    self.scale = 1.0
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_01 = buttons_panel(self)
    self.panel_02 = multi_img_scrollable(self, self.bmp_lst)
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(self.panel_01, 0, wx.EXPAND)
    sizer.Add(self.panel_02, 1, wx.EXPAND)
    self.SetSizer(sizer)
    self.Show(True)


  def _mi_list_of_wxbitmaps(self):
    bmp_obj = wxbitmap_convert(self.flex_arr)
    return bmp_obj.get_wxbitmap_lst(show_nums = self.show_nums, scale = self.scale)


  def _to_hide_nums(self):
    self.show_nums = False
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def _to_show_nums(self):
    self.show_nums = True
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def _to_re_zoom(self, rot_sn):
    if( rot_sn > 0 ):
      self.scale = self.scale * 1.05

    elif( rot_sn < 0):
      self.scale = self.scale * 0.95

    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


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
