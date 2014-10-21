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

import math

class show_3d(object):
  def __init__(self, flex_arr_in):




    if type(flex_arr_in) is list:
        print 'a list is not implemented yet'
    else:
        print 'a flex'
        app = show_3d_wx_app(redirect=False)
        app.in_lst(flex_arr_in)
        app.MainLoop()


class show_3d_wx_app(wx.App):
  def OnInit(self):
    self.frame = flex_3d_frame(None, '3D flex array viewer')
    self.panel = flex_arr_3d_outer_panel(self.frame)
    self.frame.frame_ini_img(self.panel)
    return True


  def in_lst(self, flex_lst):
    self.panel.ini_n_intro(flex_lst)
    self.SetTopWindow(self.frame)
    self.frame.Show()


class flex_3d_frame(wx.Frame):
  def __init__(self, parent, title):
    super(flex_3d_frame, self).__init__(parent, title = title,
          size = wx.DefaultSize)

  def frame_ini_img(self, in_panel):
    self.my_panel = in_panel
    self.my_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.my_sizer.Add(self.my_panel, 1, wx.EXPAND)
    self.SetSizer(self.my_sizer)


class flex_arr_3d_outer_panel(wx.Panel):
  def __init__(self, parent_frame):
    super(flex_arr_3d_outer_panel, self).__init__(parent_frame)
    self.show_nums = True


  def ini_n_intro(self, flex_arr_in):

    self.flex_arr = flex_arr_in
    self.scale = 1.0
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_01 = buttons_panel(self)
    self.panel_02 = multi_img_scrollable(self, self.bmp_lst)
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(self.panel_01, 0, wx.ALIGN_CENTRE)
    sizer.Add(self.panel_02, 1, wx.EXPAND)
    self.SetSizer(sizer)
    self.Show(True)


  def _mi_list_of_wxbitmaps(self):
    bmp_obj = wxbitmap_convert(self.flex_arr)
    return bmp_obj.get_wxbitmap_lst(show_nums = self.show_nums,
                                    scale = self.scale)


  def _to_hide_nums(self):
    self.show_nums = False
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def to_show_nums(self):
    self.show_nums = True
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def to_re_zoom(self, rot_sn):
    if( rot_sn > 0 ):
      for ntimes in range(int(math.fabs(rot_sn))):
        self.scale = self.scale * 1.05
        if( self.scale > 3.0 ):
          self.scale = 3.0
          print "Reached maximum possible Zoom"


    elif( rot_sn < 0):
      for ntimes in range(int(math.fabs(rot_sn))):
        self.scale = self.scale * 0.95
        if( self.scale < 0.2 ):
          self.scale = 0.2
          print "Reached minimum possible Zoom"


    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)

class multi_img_scrollable(scroll_pan.ScrolledPanel):
  def __init__(self, outer_panel, bmp_lst_in):
    super(multi_img_scrollable, self).__init__(outer_panel)
    self.parent_panel  = outer_panel
    self.local_bmp_lst = bmp_lst_in
    self.set_scroll_content()
    self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
    self.Bind(wx.EVT_IDLE, self.OnIdle)
    self.SetupScrolling()

    self.scroll_rot = 0


  def set_scroll_content(self):

    self.img_lst_sizer = wx.BoxSizer(wx.HORIZONTAL)
    for i, bmp_lst in enumerate(self.local_bmp_lst):
      local_bitmap = wx.StaticBitmap(self, bitmap = bmp_lst)
      slice_string = "Slice[" + str(i) + ":" + str(i + 1) + ", :, :]"
      data_txt_01 = wx.StaticText(self, -1, slice_string)
      sigle_slice_sizer = wx.BoxSizer(wx.VERTICAL)
      sigle_slice_sizer.Add(local_bitmap, wx.ALIGN_CENTRE | wx.ALL, border = 4)
      sigle_slice_sizer.Add(data_txt_01, wx.ALL, border = 4)
      self.img_lst_sizer.Add(sigle_slice_sizer,
                             flag=wx.ALIGN_CENTER | wx.ALL, border = 4)

    self.SetSizer(self.img_lst_sizer)


  def OnMouseWheel(self, event):
    sn_mov = math.copysign(1, float(event.GetWheelRotation()))
    self.scroll_rot += sn_mov


  def img_refresh(self, bmp_lst_new):
    self.local_bmp_lst = bmp_lst_new
    for child in self.GetChildren():
      child.Destroy()
    self.set_scroll_content()
    self.Layout()
    self.parent_panel.Layout()
    self.Refresh()


  def OnIdle(self, event):
    if( self.scroll_rot != 0 ):
      self.parent_panel.to_re_zoom(self.scroll_rot)
      self.scroll_rot = 0


class buttons_panel(wx.Panel):
  def __init__(self, outer_panel):
    super(buttons_panel, self).__init__(outer_panel)
    self.parent_panel  = outer_panel

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
    self.parent_panel._to_hide_nums()


  def OnShwIBut(self, event):
    self.parent_panel.to_show_nums()

