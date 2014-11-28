#!/usr/bin/env python
#
#  viewer_low_level_util.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from dials.array_family import flex
from from_flex_to_wxbitmap import wxbitmap_convert

import wx
import numpy as np
import wx.lib.scrolledpanel as scroll_pan
import wx.grid as gridlib
import math

class flex_3d_frame(wx.Frame):
  def __init__(self, parent, title):
    super(flex_3d_frame, self).__init__(parent, title = title,
          size = wx.DefaultSize)


  def frame_ini_img(self, in_upper_panel, text_data = None):
    self.img_panel = in_upper_panel

    if( text_data != None):
      self.myGrid = text_data

    self.my_sizer = wx.BoxSizer(wx.VERTICAL)
    self.my_sizer.Add(self.img_panel, proportion = 1,
                      flag =  wx.EXPAND | wx.ALL, border = 3)

    if( text_data != None):
      self.my_sizer.Add(self.myGrid, proportion = 1,
                        flag =  wx.EXPAND | wx.ALL, border = 3)

    self.my_sizer.SetMinSize((50, 20))

    self.SetSizer(self.my_sizer)



class MyGrid(gridlib.Grid):

  def __init__(self, parent_frame):
    """Constructor"""
    super(MyGrid, self).__init__(parent_frame)

  def ini_n_intro(self, table_in):
    print "Data to be visualized:"
    print "_____________________"

    col_names = []

    for col_num, col_key in enumerate(table_in[0]):
      print col_num, col_key
      if( col_key != 'shoebox' ):
        col_names.append(col_key)
    print "num of col to show =", len(col_names)

    lst_nm = range(1, 20)
    info_lst = []

    for nm in lst_nm:
      info_lst.append(table_in[nm]['miller_index'])

    print "info_lst ="
    print info_lst
    output_full_row = '''

    table[0] = {
    'imageset_id': 0,
    'xyzcal.mm': (0.0, 0.0, 0.0),
    'intensity.sum.value': 15.0,
    'xyzobs.px.variance': (0.286144578313253, 0.2676706827309237, 0.25),
    'intensity.sum.variance': 15.0,
    's1': (-0.03497287541155274, 0.604366929910932, -0.826295905178965),
    'shoebox': <dials_model_data_ext.Shoebox object at 0x7b4f820>,
    'iobs': 0,
    'rlp': (-0.035319267434222125, 0.27946863045406983, -0.5712908827798997),
    'xyzobs.mm.variance': (0.00846530120481, 0.007918769477911, 1.7134729863002e-06),
    'miller_index': (36, -4, -17),
    'flags': 0,
    'bbox': (1158, 1161, 68, 70, 0, 1),
    'xyzobs.mm.value': (199.43804512406385, 11.977194484339728, 1.4324789835743459),
    'xyzobs.px.value': (1159.5, 69.23333333333333, 0.5),
    'id': 0,
    'panel': 0
    }

    '''


    self.CreateGrid(len(info_lst), len(col_names))
    #for nm in lst_nm:
    for nm, data in enumerate(info_lst):
      print "nm =", nm
      #self.SetCellValue(nm - 1, 3, str(info_lst[nm - 1]))
      for col_pos, cel_val in enumerate( table_in[nm].iteritems() ):

        #self.SetCellValue(nm, 3, str(data))
        self.SetCellValue(nm, col_pos, str(cel_val))

    self.EnableEditing(False)
    '''self.SetCellSize(row(pos), col(pos),
                      size(n of grid rows), size(n of grid cols))'''
    #self.SetCellSize(5, 3, 2, 50)
    self.AutoSizeColumns(True)
    self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)

  def OnCellLeftClick(self, evt):
    print "OnCellLeftClick: (%d,%d) %s\n" % (evt.GetRow(),
                                             evt.GetCol(),
                                             evt.GetPosition())
    evt.Skip()

class flex_arr_img_panel(wx.Panel):
  def __init__(self, parent_frame):
    super(flex_arr_img_panel, self).__init__(parent_frame)
    self.show_nums = True
    self.show_mask = True

  def ini_n_intro(self, flex_arr_one, flex_arr_two = None):
    self.first_lst_in, self.segn_lst_in = flex_arr_one, flex_arr_two
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
      if(self.show_mask == True):
        self.lst_bmp_obj = wxbitmap_convert(self.first_lst_in, self.segn_lst_in)
      else:
        self.lst_bmp_obj = wxbitmap_convert(self.first_lst_in, None)
      return self.lst_bmp_obj.get_wxbitmap_lst(show_nums = self.show_nums,
                                      scale = self.scale)

    else:
      return self.lst_bmp_obj.scaling(scale = self.scale)


  def to_hide_nums(self):
    self.show_nums = False
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def to_show_nums(self):
    self.show_nums = True
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def to_show_mask(self):
    self.show_mask = True
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def to_hide_mask(self):
    self.show_mask = False
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)


  def to_re_zoom(self, rot_sn):
    if( rot_sn > 0 ):
      for ntimes in range(int(math.fabs(rot_sn))):
        self.scale = self.scale * 1.05
        if( self.scale > 3.0 ):
          #Maximum possible zoom reached
          self.scale = 3.0

    elif( rot_sn < 0):
      for ntimes in range(int(math.fabs(rot_sn))):
        self.scale = self.scale * 0.95
        if( self.scale < 0.2 ):
          #Minimum possible zoom reached
          self.scale = 0.2

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
    self.SetBackgroundColour(wx.Colour(200, 200, 200))
    aprox_len_pix = len(self.lst_2d_bmp) * 10
    self.SetScrollbars(1, 1, aprox_len_pix * 10, aprox_len_pix * 10)


  def set_scroll_content(self):

    img_lst_vert_sizer = wx.BoxSizer(wx.VERTICAL)

    for lst_1d in self.lst_2d_bmp:
      img_lst_hor_sizer = wx.BoxSizer(wx.HORIZONTAL)

      for i, bmp_lst in enumerate(lst_1d):
        local_bitmap = wx.StaticBitmap(self, bitmap = bmp_lst)
        slice_string = "Slice[" + str(i) + ":" + str(i + 1) + ", :, :]"
        slice_sub_info_txt = wx.StaticText(self, -1, slice_string)
        sigle_slice_sizer = wx.BoxSizer(wx.VERTICAL)
        sigle_slice_sizer.Add(local_bitmap, proportion = 0,
                              flag = wx.ALIGN_CENTRE | wx.ALL, border = 2)
        sigle_slice_sizer.Add(slice_sub_info_txt, proportion = 0,
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

    Show_Its_CheckBox = wx.CheckBox(self, -1, "Show I nums")
    Show_Its_CheckBox.SetValue(True)
    Show_Its_CheckBox.Bind(wx.EVT_CHECKBOX, self.OnItsCheckbox)

    if( self.parent_panel.segn_lst_in != None ):
      Show_Msk_CheckBox = wx.CheckBox(self, -1, "Show Mask")
      Show_Msk_CheckBox.SetValue(True)
      Show_Msk_CheckBox.Bind(wx.EVT_CHECKBOX, self.OnMskCheckbox)

    self.my_sizer = wx.BoxSizer(wx.VERTICAL)
    self.my_sizer.Add(Show_Its_CheckBox, proportion = 0,
                      flag = wx.ALIGN_TOP, border = 5)

    if( self.parent_panel.segn_lst_in != None ):
      self.my_sizer.Add(Show_Msk_CheckBox, proportion = 0,
                        flag = wx.ALIGN_TOP, border = 5)

    self.my_sizer.SetMinSize((50, 20))
    self.SetSizer(self.my_sizer)


  def OnItsCheckbox(self, event):
    debugg_screen_log = '''
    print "OnItsCheckbox"
    print "event.IsChecked() =", event.IsChecked()
    '''
    if(event.IsChecked() == True):
      self.parent_panel.to_show_nums()
    else:
      self.parent_panel.to_hide_nums()


  def OnMskCheckbox(self, event):
    debugg_screen_log = '''
    print "OnMskCheckbox"
    print "event.IsChecked() =", event.IsChecked()
    '''
    if(event.IsChecked() == True):
      self.parent_panel.to_show_mask()
    else:
      self.parent_panel.to_hide_mask()

