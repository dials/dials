#!/usr/bin/env python
#
#  multi_3D_slice_viewer.py
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

#self.CalcScrolledPosition((old_x,old_y))

class show_3d(object):
  def __init__(self, flex_arr_in):
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


  def _mi_list_of_wxbitmaps(self, re_scaling = False):
    if(re_scaling == False):
      self.lst_bmp_obj = wxbitmap_convert(self.flex_arr)
      return self.lst_bmp_obj.get_wxbitmap_lst(show_nums = self.show_nums,
                                      scale = self.scale)

    else:
      return self.lst_bmp_obj.scaling(scale = self.scale)


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
          print "Maximum possible zoom reached "

    elif( rot_sn < 0):
      for ntimes in range(int(math.fabs(rot_sn))):
        self.scale = self.scale * 0.95
        if( self.scale < 0.2 ):
          self.scale = 0.2
          print "Minimum possible zoom reached"


    self.bmp_lst = self._mi_list_of_wxbitmaps(re_scaling = True)
    self.panel_02.img_refresh(self.bmp_lst)

class multi_img_scrollable(scroll_pan.ScrolledPanel):
  def __init__(self, outer_panel, bmp_lst_in):
    super(multi_img_scrollable, self).__init__(outer_panel)
    self.parent_panel  = outer_panel
    self.lst_2d_bmp = bmp_lst_in
    self.set_scroll_content()
    self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
    self.Bind(wx.EVT_IDLE, self.OnIdle)
    self.SetupScrolling()
    self.scroll_rot = 0
    self.SetBackgroundColour(wx.Colour(200,200,200))

    aprox_len_pix = len(self.lst_2d_bmp) * 10
    #print "aprox_len_pix =", aprox_len_pix
    self.SetScrollbars(1, 1, aprox_len_pix * 10, aprox_len_pix * 10)

  def set_scroll_content(self):

    img_lst_vert_sizer = wx.BoxSizer(wx.VERTICAL)

    for lst_1d in self.lst_2d_bmp:
      img_lst_hor_sizer = wx.BoxSizer(wx.HORIZONTAL)

      for i, bmp_lst in enumerate(lst_1d):
        local_bitmap = wx.StaticBitmap(self, bitmap = bmp_lst)
        slice_string = "Slice[" + str(i) + ":" + str(i + 1) + ", :, :]"
        data_txt_01 = wx.StaticText(self, -1, slice_string)
        sigle_slice_sizer = wx.BoxSizer(wx.VERTICAL)
        sigle_slice_sizer.Add(local_bitmap, proportion = 0,
                              flag = wx.ALIGN_CENTRE | wx.ALL, border = 2)
        sigle_slice_sizer.Add(data_txt_01, proportion = 0,
                              flag = wx.ALIGN_CENTRE | wx.ALL, border = 2)
        img_lst_hor_sizer.Add(sigle_slice_sizer, proportion = 0,
                              flag = wx.ALIGN_CENTER | wx.ALL, border = 2)

      img_lst_vert_sizer.Add(img_lst_hor_sizer, proportion = 0,
                             flag = wx.ALIGN_CENTER | wx.TOP, border = 6)

    self.SetSizer(img_lst_vert_sizer)





  def OnMouseWheel(self, event):

    #saving amount of scroll steps to do
    sn_mov = math.copysign(1, float(event.GetWheelRotation()))
    self.scroll_rot += sn_mov

    #saving relative position of scrolling area to keep after scrolling
    v_size_x, v_size_y = self.GetVirtualSize()

    self.x_to_keep, self.y_to_keep = self.GetViewStart()
    self.x_to_keep = float(self.x_to_keep) / float(v_size_x)
    self.y_to_keep = float(self.y_to_keep) / float(v_size_y)

    debugg_screen_log = '''
    print "self.GetVirtualSize() =", self.GetVirtualSize()
    print "self.GetViewStart() =", self.GetViewStart()
    print "self.GetScrollPixelsPerUnit() =", self.GetScrollPixelsPerUnit()
    print "self.x_to_keep =", self.x_to_keep
    print "self.y_to_keep =", self.y_to_keep
    '''

  def img_refresh(self, bmp_lst_new):
    self.lst_2d_bmp = bmp_lst_new
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
      v_size_x, v_size_y = self.GetVirtualSize()
      self.Scroll(self.x_to_keep * v_size_x, self.y_to_keep * v_size_y)


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

