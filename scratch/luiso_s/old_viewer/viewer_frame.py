#
#  DIALS viewer_frame
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

import wx
from wx.lib.buttons import GenBitmapTextButton
from dials.scratch.luiso_s.old_viewer.viewer_utilities import np_to_bmp

from dials.scratch.luiso_s.old_viewer.reflection_data_navigator import table_s_navigator

class ReflectionFrame(wx.Frame):
  def __init__(self, *args, **kwargs):
    wx.Frame.__init__(self, *args, **kwargs)

    #self.MaxImageSizeX = 280
    #self.MaxImageSizeY = 150
    self.MaxImageSizeX = 320
    self.MaxImageSizeY = 260

    lbl_ref = "Reflection"

    img_tmp = wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD, wx.ART_OTHER, (16, 16))
    btn_nxt_refl = wx.lib.buttons.GenBitmapTextButton(self
                   , label = lbl_ref, bitmap = img_tmp)
    img_tmp = wx.ArtProvider.GetBitmap(wx.ART_GO_BACK, wx.ART_OTHER, (16, 16))
    btn_prv_refl = wx.lib.buttons.GenBitmapTextButton(self
                   , label = lbl_ref, bitmap = img_tmp)
    btn_read = wx.Button(self, -1, "Reflection #")

    radio1 = wx.RadioButton(self, -1, "data, background, mask"
                           , style = wx.RB_GROUP)
    radio2 = wx.RadioButton(self, -1, "3 layers of data")
    radio3 = wx.RadioButton(self, -1, "3 layers of background")
    radio4 = wx.RadioButton(self, -1, "3 layers of mask")

    img_tmp = wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD, wx.ART_OTHER, (16, 16))
    btn_nxt_slice = wx.lib.buttons.GenBitmapTextButton(self
                   , label = 'slice', bitmap = img_tmp)

    img_tmp = wx.ArtProvider.GetBitmap(wx.ART_GO_BACK, wx.ART_OTHER, (16, 16))
    btn_prv_slice = wx.lib.buttons.GenBitmapTextButton(self
                   , label = 'slice', bitmap = img_tmp)

    self.data_txt_01 = wx.StaticText(self, -1, "(data_txt)", size = (800, 16))
    self.data_txt_02 = wx.StaticText(self, -1, "(data_txt)", size = (800, 16))
    self.data_txt_03 = wx.StaticText(self, -1, "(data_txt)", size = (800, 16))

    self.Image = [None, None, None]

    for indx in range(len(self.Image)):
      self.Image[indx] = wx.StaticBitmap(self, bitmap = wx.EmptyBitmap(
                                         self.MaxImageSizeX, self.MaxImageSizeY))

    self.text01 = wx.TextCtrl(self, -1, "enter number", size = (100, -1))

    u_box = wx.BoxSizer(wx.VERTICAL)
    div_h_box = wx.BoxSizer(wx.HORIZONTAL)
    div_h_box.Add(btn_prv_refl, 0, wx.CENTER | wx.ALL,3)
    div_h_box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,3)
    u_box.Add(div_h_box)

    div_h_box = wx.BoxSizer(wx.HORIZONTAL)
    div_h_box.Add(btn_read, 0, wx.CENTER | wx.ALL,3)
    div_h_box.Add(self.text01, 0, wx.CENTER | wx.ALL,3)
    u_box.Add(div_h_box)

    u_box.Add(radio1, 0, wx.LEFT | wx.ALL,3)
    u_box.Add(radio2, 0, wx.LEFT | wx.ALL,3)
    u_box.Add(radio3, 0, wx.LEFT | wx.ALL,3)
    u_box.Add(radio4, 0, wx.LEFT | wx.ALL,3)

    div_h_box = wx.BoxSizer(wx.HORIZONTAL)
    div_h_box.Add(btn_prv_slice, 0, wx.CENTER | wx.ALL,3)
    div_h_box.Add(btn_nxt_slice, 0, wx.CENTER | wx.ALL,3)
    u_box.Add(div_h_box)

    h_box = wx.BoxSizer(wx.HORIZONTAL)

    h_box.Add(u_box)

    div_h_box = wx.BoxSizer(wx.HORIZONTAL)
    div_h_box.Add(self.Image[0], 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 5)

    div_h_box.Add(self.Image[1], 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 5)

    div_h_box.Add(self.Image[2], 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 5)
    div_u_box = wx.BoxSizer(wx.VERTICAL)
    div_u_box.Add(div_h_box)
    div_u_box.Add(self.data_txt_01, 0, wx.CENTER | wx.ALL,3)
    div_u_box.Add(self.data_txt_02, 0, wx.CENTER | wx.ALL,3)
    div_u_box.Add(self.data_txt_03, 0, wx.CENTER | wx.ALL,3)
    h_box.Add(div_u_box)

    self.frame_scale = 0.4
    self.opt = 0

    self.sizing_counter = 0
    self.SetSizerAndFit(h_box)

    btn_nxt_refl.Bind(wx.EVT_BUTTON, self._DisplayNext_refl)
    btn_prv_refl.Bind(wx.EVT_BUTTON, self._DisplayPrev_refl)
    btn_nxt_slice.Bind(wx.EVT_BUTTON, self._DisplayNext_slice)
    btn_prv_slice.Bind(wx.EVT_BUTTON, self._DisplayPrev_slice)

    self.Bind(wx.EVT_RADIOBUTTON, self._OnRadio1, radio1)
    self.Bind(wx.EVT_RADIOBUTTON, self._OnRadio2, radio2)
    self.Bind(wx.EVT_RADIOBUTTON, self._OnRadio3, radio3)
    self.Bind(wx.EVT_RADIOBUTTON, self._OnRadio4, radio4)

    self.Bind(wx.EVT_BUTTON, self._read_num_from_txt, btn_read)

    self.Bind(wx.EVT_SIZE, self._OnSize)

    self.bmp = np_to_bmp()
    self.arr_img = [None, None, None]
    self.arr_title = [None, None, None]
    self.My_Img = [None, None, None]
    self.wx_Img = [None, None, None]
    wx.EVT_CLOSE(self, self._On_Close_Window)

  def tabl_to_frame(self, loc_tabl):
    self.tabl = table_s_navigator(loc_tabl)
    self._DisplayPrev_refl()
  def _DisplayNext_refl(self, event = None):
    self.tabl.next_Reflection()
    self._My_Update()
  def _DisplayPrev_refl(self, event = None):
    self.tabl.Previous_Reflection()
    self._My_Update()
  def _DisplayNext_slice(self, event = None):
    self.tabl.next_slice()
    self._My_Update()
  def _DisplayPrev_slice(self, event = None):
    self.tabl.Previous_slice()
    self._My_Update()

  def _read_num_from_txt(self, event = None):
    # TODO check if the user entered a number
    a = int(self.text01.GetValue())
    self.tabl.Jump_to_Reflection(a)
    self._My_Update()

  def _OnRadio1(self, event = None):
    self.opt = 0
    # opt = 0 means (data, background, mask)
    self._My_Update(request_new_data = True)

  def _OnRadio2(self, event = None):
    self.opt = 1
    # self.opt = 1 means (3 layers of data)
    self._My_Update(request_new_data = True)

  def _OnRadio3(self, event = None):
    self.opt = 2
    # self.opt = 2 means (3 layers of background)
    self._My_Update(request_new_data = True)

  def _OnRadio4(self, event = None):
    self.opt = 3
    # self.opt = 3 means (3 layers of mask)
    self._My_Update(request_new_data = True)

  def _OnSize(self, event = None):
    if( self.sizing_counter > 5 ):
      siz_data = self.GetSize()
      #print "New size of window =", siz_data
      optm_aspec_ratio = 4.21
      #print "siz_data = ", siz_data[0], siz_data[1]
      #print "aspect ratio = ", float(siz_data[0])/ float(siz_data[1])
      aspec_ratio = float(siz_data[0])/ float(siz_data[1])

      if(aspec_ratio > optm_aspec_ratio):
        #print "use float(siz_data[1] (Height) to calculate new size"
        self.frame_scale = float(siz_data[1]) * 0.5 / 295.0
        #(1227, 295)
      else:
        #print "use float(siz_data[0] (with) to calculate new size"
        self.frame_scale = float(siz_data[0]) * 0.5 / 1227.0
        #(1227, 295)
      self._My_Update(request_new_data = False)
      #print "resizing"
      #print "aspec_ratio =", aspec_ratio
    else:
      self.sizing_counter += 1

  def _My_Update(self, request_new_data = True):

    if( request_new_data == True ):
      self.arr_img, self.arr_title = self.tabl(opt = self.opt)
      ref_max = self.tabl.Get_Max(self.opt)
      #box_lmt = self.tabl.Get_bbox()
      self.xyz_px = self.tabl.Get_xyz()
      self.bbox_px = self.tabl.Get_bbox()
      self.hkl_data = self.tabl.Get_hkl()
      for indx in range(len(self.arr_img)):
        if( (self.opt == 0 and indx == len(self.arr_img) - 1) or self.opt == 3):
          Imax = 10
          mask_colour = True
        else:
          Imax = ref_max
          mask_colour = False
        #print "calling np_to_bmp.__init__(...)"
        self.wx_Img[indx] = self.bmp(np_img_2d = self.arr_img[indx],
                                     Intst_max = Imax, ofst = self.bbox_px,
                                     xyz = self.xyz_px, title = self.arr_title[indx],
                                     mask_colour = mask_colour)
      self.text01.SetValue(str(self.tabl.Get_ref_num()))

    for indx in range(len(self.arr_img)):
      if( self.arr_img[indx] == None ):
        empty_flag = True
      else:
        empty_flag = False
      self.My_Img[indx] = self.bmp.from_wx_image_to_wx_bitmap(self.wx_Img[indx]
                        , self.frame_scale, empty = empty_flag)
      self.Image[indx].SetBitmap(self.My_Img[indx])
      if( self.xyz_px != None ):
        my_str = ( " Centroid Position = ( "
                 + str(self.xyz_px[0]) + ", "
                 + str(self.xyz_px[1]) + ", "
                 + str(self.xyz_px[2]) + ") ")
        self.data_txt_01.SetLabel(my_str)
      else:
        self.data_txt_01.SetLabel(" No (x, y, z) Data")

      if( self.bbox_px != None ):
        my_str_01 = ( " Bbox = ( "
                 + str(self.bbox_px[0]) + ", "
                 + str(self.bbox_px[1]) + ", "
                 + str(self.bbox_px[2]) + ", "
                 + str(self.bbox_px[3]) + ", "
                 + str(self.bbox_px[4]) + ", "
                 + str(self.bbox_px[5]) + ") ")
        self.data_txt_02.SetLabel(my_str_01)
      else:
        self.data_txt_02.SetLabel(" No Bbox Data")

      try:

        hkl_string = 'HKL:    ' + str(self.hkl_data[0]) + \
                            '    ' + str(self.hkl_data[1]) + \
                            '    ' + str(self.hkl_data[2]) + ' '
      except:
        hkl_string = ' No HKL data'
      self.data_txt_03.SetLabel(hkl_string)

    self.Layout()
    self.Refresh()

  def _On_Close_Window(self, event):
    self.Destroy()
