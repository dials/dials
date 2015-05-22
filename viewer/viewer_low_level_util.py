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
from __future__ import division
from dials.array_family import flex
from from_flex_to_wxbitmap import wxbitmap_convert

import wx
import wx.lib.scrolledpanel as scroll_pan
import wx.grid as gridlib
import math


class grid_frame(wx.Frame):
  def __init__(self, parent, title):
    super(grid_frame, self).__init__(parent, title = title,
          size = wx.DefaultSize)

  def frame_ini_img(self, in_upper_panel, text_data = None):
    self.img_panel = in_upper_panel



    if( text_data != None ):
      self.myGrid = text_data

    self.my_sizer = wx.BoxSizer(wx.VERTICAL)

    if( text_data != None ):
      self.my_sizer.Add(self.myGrid, proportion = 3,
                        flag =  wx.EXPAND | wx.ALL, border = 3)



    self.my_sizer.SetMinSize((900, 600))
    self.SetSizer(self.my_sizer)
    self.my_sizer.Fit(self)
    self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

  def OnCloseWindow(self, event):
    wx.GetApp().Exit()



class flex_3d_frame(wx.Frame):
  def __init__(self, parent, title):
    super(flex_3d_frame, self).__init__(parent, title = title,
          size = wx.DefaultSize)
    self.table_exist = False

  def frame_ini_img(self, in_upper_panel, text_data = None):
    self.img_panel = in_upper_panel
    self.myGrid = None

    if( text_data != None ):
      self.myGrid = text_data


    self.my_sizer = wx.BoxSizer(wx.VERTICAL)
    self.my_sizer.Add(self.img_panel, proportion = 2,
                      flag =  wx.EXPAND | wx.ALL, border = 3)

    if( text_data != None ):
      self.my_sizer.Add(self.myGrid, proportion = 3,
                        flag =  wx.EXPAND | wx.ALL, border = 3)


    self.my_sizer.SetMinSize((400, 200))
    self.SetSizer(self.my_sizer)
    self.my_sizer.Fit(self)
    self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

  def OnCloseWindow(self, event):
    if( self.table_exist == False ):
      #print "from flex_3d_frame self.myGrid = None"
      event.Skip()
    else:
      #print "from flex_3d_frame self.myGrid .NEQ. None"
      wx.GetApp().Exit()


class TupTable(gridlib.PyGridTableBase):
  def __init__(self, data, rowLabels=None, colLabels=None):
    gridlib.PyGridTableBase.__init__(self)
    self.data = data
    self.rowLabels = rowLabels
    self.colLabels = colLabels

  def GetNumberRows(self):
    return len(self.data)

  def GetNumberCols(self):
    return len(self.data[0])

  def GetColLabelValue(self, col):
    if self.colLabels:
      return self.colLabels[col]

  def GetRowLabelValue(self, row):
    if self.rowLabels:
      return self.rowLabels[row]

  def IsEmptyCell(self, row, col):
    return False

  def GetValue(self, row, col):
    return self.data[row][col]

  def SetValue(self, row, col, value):
    pass


class MyGrid(gridlib.Grid):

  def __init__(self, parent_frame):

    self.parent_fr = parent_frame
    super(MyGrid, self).__init__(parent_frame)

  def ini_n_intro(self, table_in):

    self.lst_keys = []
    self.data = []
    self.sorted_flags = []

    for key in table_in.keys():
      if(key != "shoebox"):
        col_label = str(key)
        col_content = map(str, table_in[key])

        str_len = bigger_size(col_label, col_content)

        if(str_len > len(col_label)):
          lng_dif = int(float(str_len - len(col_label)) * 1.6)
          add_str = ' ' * lng_dif
          col_label = add_str + col_label + add_str
        else:
          col_label = "   " + col_label + "   "

        self.lst_keys.append(col_label)
        self.data.append(col_content)
        self.sorted_flags.append(True)

    self.lst_keys.append("lst pos")
    self.data.append(range(len(table_in)))
    self.sorted_flags.append(True)

    self.last_col_num = len(self.lst_keys) - 1

    self.data = tuple(zip(*self.data))
    self.EnableEditing(True)
    self.EnableDragGridSize(True)

    self.set_my_table(self.last_col_num)

    self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)
    self.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.OnLabelLeftClick)

  def OnLabelLeftClick(self, evt):

    #timing_for_debugging = '''
    import time
    time1 = time.time()
    #'''

    if(evt.GetCol() == -1):
      self.repaint_img(evt.GetRow())
    else:
      self.set_my_table(evt.GetCol())

    evt.Skip()

    #timing_for_debugging = '''
    time2 = time.time()
    timedif = time2 - time1
    print "timedif =", timedif
    #'''

  def set_my_table(self, col_to_sort):

    self.sorted_flags[col_to_sort] = not self.sorted_flags[col_to_sort]



    try:
      tupldata = sorted(self.data, key=lambda x: int(x[col_to_sort]),
                        reverse = self.sorted_flags[col_to_sort])
      print "using key=lambda x: int(x[col_to_sort])"

    except:
      try:
        tupldata = sorted(self.data, key=lambda x: float(x[col_to_sort]),
                          reverse = self.sorted_flags[col_to_sort])
        print "using key=lambda x: float(x[col_to_sort])"
      except:
        tupldata = sorted(self.data, key=lambda x:
                          tuple(eval(str(x[col_to_sort]))),
                          reverse = self.sorted_flags[col_to_sort])
        print "using key=lambda x: tuple(x[col_to_sort])"

    colLabels = tuple(self.lst_keys)
    rowLabels = tuple(range(len(tupldata)))

    tableBase = TupTable(tupldata, rowLabels, colLabels)

    #self.AutoSizeColumns(False)
    self.SetTable(tableBase)
    self.Refresh()
    #self.AutoSizeColumn(1)
    for i in range( len(self.lst_keys) ):
      self.AutoSizeColLabelSize(i)

    #self.AutoSizeColumns(True)

  def OnCellLeftClick(self, evt):
    self.repaint_img(evt.GetRow())
    evt.Skip()

  def repaint_img(self, clikd_row):
    new_row = int(self.GetCellValue(clikd_row, self.last_col_num))
    self.parent_fr.img_panel.update_img_w_row_pos(new_row)

class flex_arr_img_panel(wx.Panel):
  def __init__(self, parent_frame):
    super(flex_arr_img_panel, self).__init__(parent_frame)
    self.show_nums = True
    self.show_mask = True
    self.row_pos = 0
    self.Pframe = parent_frame

  def ini_n_intro(self, data_in_one, data_in_two = None):

    self.scale = 1.0

    if( isinstance(data_in_one, flex.reflection_table) ):
      self.table = data_in_one
      self.assign_row_pos()

    else:
      self.first_lst_in, self.segn_lst_in = data_in_one, data_in_two

    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_01 = buttons_panel(self)
    self.panel_02 = multi_img_scrollable(self, self.bmp_lst)
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(self.panel_01, 0, wx.EXPAND)
    sizer.Add(self.panel_02, 1, wx.EXPAND)
    self.SetSizer(sizer)
    self.Show(True)

  def assign_row_pos(self):

    try:
      self.first_lst_in, \
      self.segn_lst_in = \
                              self.table[self.row_pos]['shoebox'].data, \
                              self.table[self.row_pos]['shoebox'].mask
    except:
      self.first_lst_in, self.segn_lst_in = None, None

  def _mi_list_of_wxbitmaps(self, re_scaling = False):

    if( re_scaling == False ):
      if( self.show_mask == True ):
        self.lst_bmp_obj = wxbitmap_convert(self.first_lst_in, self.segn_lst_in)
      else:
        self.lst_bmp_obj = wxbitmap_convert(self.first_lst_in, None)
      return self.lst_bmp_obj.get_wxbitmap_lst(show_nums = self.show_nums,
                                      scale = self.scale)

    else:
      return self.lst_bmp_obj.scaling(scale = self.scale)

  def update_img_w_row_pos(self, num):
    self.row_pos = num
    self.assign_row_pos()
    self.bmp_lst = self._mi_list_of_wxbitmaps()
    self.panel_02.img_refresh(self.bmp_lst)

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
  def __init__(self, outer_panel, i_bmp_in):
    super(multi_img_scrollable, self).__init__(outer_panel)
    self.parent_panel  = outer_panel
    self.lst_2d_bmp = i_bmp_in
    self.SetBackgroundColour(wx.Colour(200, 200, 200))

    self.mainSizer = wx.BoxSizer(wx.VERTICAL)

    self.n_img = 0
    self.img_lst_v_sizer = wx.BoxSizer(wx.VERTICAL)

    self.set_scroll_content()

    self.mainSizer.Add(self.img_lst_v_sizer, 0, wx.CENTER|wx.ALL, 10)
    self.SetSizer(self.mainSizer)

    self.SetupScrolling()
    self.SetScrollRate(1, 1)
    self.Mouse_Pos_x = -1
    self.Mouse_Pos_y = -1
    self.old_Mouse_Pos_x = -1
    self.old_Mouse_Pos_y = -1
    self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
    self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftButDown)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftButUp)
    self.Bind(wx.EVT_IDLE, self.OnIdle)
    self.scroll_rot = 0

    not_working = '''

    for child in self.GetChildren():
      child.Bind(wx.EVT_MOUSE, self.OnLeftButDown)

    for lst_1d in self.lst_2d_bmp:
      for i_bmp in lst_1d:
        print "i_bmp =", i_bmp
        #i_bmp.Bind(wx.EVT_LEFT_DOWN, self.OnLeftButDown)
    '''

  def OnLeftButDown(self, event):
    #print "OnLeftButDown(self)"
    self.Bdwn = True
    self.old_Mouse_Pos_x, self.old_Mouse_Pos_y = event.GetPosition()

  def OnLeftButUp(self, event):
    self.Bdwn = False
    Mouse_Pos_x, Mouse_Pos_y = event.GetPosition()
    #print "Mouse_Pos_x, Mouse_Pos_y =", Mouse_Pos_x, Mouse_Pos_y

  def set_scroll_content(self):

    self.img_lst_v_sizer.Clear(True)

    for lst_1d in self.lst_2d_bmp:
      img_lst_hor_sizer = wx.BoxSizer(wx.HORIZONTAL)

      for i, i_bmp in enumerate(lst_1d):
        local_bitmap = wx.StaticBitmap(self, bitmap = i_bmp)
        slice_string = "Slice[" + str(i) + ":" + str(i + 1) + ", :, :]"
        slice_sub_info_txt = wx.StaticText(self, -1, slice_string)

        sigle_slice_sizer = wx.BoxSizer(wx.VERTICAL)
        sigle_slice_sizer.Clear(True)

        sigle_slice_sizer.Add(local_bitmap, proportion = 0,
                              flag = wx.ALIGN_CENTRE | wx.ALL, border = 2)
        sigle_slice_sizer.Add(slice_sub_info_txt, proportion = 0,
                              flag = wx.ALIGN_CENTRE | wx.ALL, border = 2)
        img_lst_hor_sizer.Add(sigle_slice_sizer, proportion = 0,
                              flag = wx.ALIGN_CENTER | wx.ALL, border = 2)

      self.img_lst_v_sizer.Add(img_lst_hor_sizer, 0, wx.CENTER|wx.ALL, 10)

      self.n_img += 1

    self.parent_panel.Pframe.Layout()
    self.SetScrollRate(1, 1)
  def OnMouseMotion(self, event):

    self.Mouse_Pos_x, self.Mouse_Pos_y = event.GetPosition()
    #print "OnMouseMotion(self, event)"



  def OnMouseWheel(self, event):

    #Getting amount of scroll steps to do
    sn_mov = math.copysign(1, float(event.GetWheelRotation()))
    self.scroll_rot += sn_mov

    #Getting relative position of scrolling area to keep after scrolling
    v_size_x, v_size_y = self.GetVirtualSize()
    self.Mouse_Pos_x, self.Mouse_Pos_y = event.GetPosition()

    View_start_x, View_start_y = self.GetViewStart()

    No_longer_needed = '''
    x_scroll_increment, y_scroll_increment = self.GetScrollPixelsPerUnit()
    print "x_scroll_increment, y_scroll_increment", x_scroll_increment, y_scroll_increment
    View_start_x = View_start_x * x_scroll_increment
    View_start_y = View_start_y * y_scroll_increment
    '''

    self.x_uni = float( View_start_x + self.Mouse_Pos_x ) / float(v_size_x)
    self.y_uni = float( View_start_y + self.Mouse_Pos_y ) / float(v_size_y)


  def img_refresh(self, i_bmp_new):

    self.lst_2d_bmp = i_bmp_new
    self.set_scroll_content()


  def OnIdle(self, event):
    if( self.scroll_rot != 0 ):

      self.SetScrollRate(1, 1)

      self.parent_panel.to_re_zoom(self.scroll_rot)
      self.scroll_rot = 0
      v_size_x, v_size_y = self.GetVirtualSize()
      new_scroll_pos_x = float(self.x_uni * v_size_x - self.Mouse_Pos_x)
      new_scroll_pos_y = float(self.y_uni * v_size_y - self.Mouse_Pos_y)

      No_longer_needed = '''
      x_scroll_increment, y_scroll_increment = self.GetScrollPixelsPerUnit()
      print "x_scroll_increment, y_scroll_increment", x_scroll_increment, y_scroll_increment
      new_scroll_pos_x = new_scroll_pos_x  / x_scroll_increment + 0.5
      new_scroll_pos_y = new_scroll_pos_y  / y_scroll_increment + 0.5
      '''
      self.Scroll(new_scroll_pos_x, new_scroll_pos_y)

    else:
      if(self.old_Mouse_Pos_x != -1 and self.old_Mouse_Pos_y != -1):
        #print "need to move IMG"
        pass


      self.old_Mouse_Pos_x = self.Mouse_Pos_x
      self.old_Mouse_Pos_y = self.Mouse_Pos_y



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

    self.my_sizer.SetMinSize((50, 300))
    self.SetSizer(self.my_sizer)

  def OnItsCheckbox(self, event):

    if(event.IsChecked() == True):
      self.parent_panel.to_show_nums()
    else:
      self.parent_panel.to_hide_nums()

  def OnMskCheckbox(self, event):

    if(event.IsChecked() == True):
      self.parent_panel.to_show_mask()
    else:
      self.parent_panel.to_hide_mask()

def bigger_size(str_label, lst_col):
  lng_label_ini = len(str_label)
  lng_lst_col = len(lst_col)
  lng_final = lng_label_ini

  if( lng_lst_col < 30 ):
    pos_lst = range(lng_lst_col)

  else:
    pos_lst = range(15)
    pos_lst += range(lng_lst_col - 15, lng_lst_col )

  lng_cel_zero = 0
  for pos in pos_lst:
    lng_cel_pos = len(str(lst_col[pos]))
    if( lng_cel_pos > lng_cel_zero ):
      lng_cel_zero = lng_cel_pos

    if(lng_cel_zero > lng_label_ini):
      lng_final = lng_cel_zero

  return lng_final
