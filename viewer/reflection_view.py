#
#  DIALS viewer
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

import wx
from dials.viewer.viewer_utilities \
     import from_wx_image_to_wx_bitmap, np_to_bmp

from dials.viewer.reflection_data_navigator import table_s_navigator

class ReflectionFrame(wx.Frame):
  def __init__(self, *args, **kwargs):
    wx.Frame.__init__(self, *args, **kwargs)


    self.MaxImageSizeX = 280
    self.MaxImageSizeY = 150

    btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
    btn_prv_refl = wx.Button(self, -1, "Previous Reflection")
    btn_nxt_slice = wx.Button(self, -1, "Next slice ")
    btn_prv_slice = wx.Button(self, -1, "Previous slice")

    radio1 = wx.RadioButton(self, -1, "data, background, mask", style = wx.RB_GROUP)
    radio2 = wx.RadioButton(self, -1, "3 layers of data")
    radio3 = wx.RadioButton(self, -1, "3 layers of background")
    radio4 = wx.RadioButton(self, -1, "3 layers of mask")

    self.Image_01 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_02 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_03 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))

    h_box = wx.BoxSizer(wx.HORIZONTAL)
    u_box = wx.BoxSizer(wx.VERTICAL)

    u_box.Add(btn_prv_refl, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(radio1, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(radio2, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(radio3, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(radio4, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(btn_nxt_slice, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(btn_prv_slice, 0, wx.CENTER | wx.ALL,5)
    h_box.Add(u_box)

    h_box.Add(self.Image_01, 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 7)

    h_box.Add(self.Image_02, 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 7)

    h_box.Add(self.Image_03, 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 7)

    self.frame_scale = 0.5
    self.opt = 0

    self.sizing_counter = 0
    self.SetSizerAndFit(h_box)

    btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext_refl)
    btn_prv_refl.Bind(wx.EVT_BUTTON, self.DisplayPrev_refl)
    btn_nxt_slice.Bind(wx.EVT_BUTTON, self.DisplayNext_slice)
    btn_prv_slice.Bind(wx.EVT_BUTTON, self.DisplayPrev_slice)

    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio1, radio1)
    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio2, radio2)
    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio3, radio3)
    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio4, radio4)

    self.Bind(wx.EVT_SIZE, self.OnSize)

    self.bmp = np_to_bmp()

    wx.EVT_CLOSE(self, self.On_Close_Window)

  def tabl_to_frame(self, loc_tabl):
    self.tabl = table_s_navigator(loc_tabl)
    self.DisplayPrev_refl()
  def DisplayNext_refl(self, event = None):
    self.tabl.next_Reflection()
    self.My_Update()
  def DisplayPrev_refl(self, event = None):
    self.tabl.Previous_Reflection()
    self.My_Update()
  def DisplayNext_slice(self, event = None):
    self.tabl.next_slice()
    self.My_Update()
  def DisplayPrev_slice(self, event = None):
    self.tabl.Previous_slice()
    self.My_Update()


  def OnRadio1(self, event = None):
    self.opt = 0
    self.My_Update(request_new_data = True)

  def OnRadio2(self, event = None):
    self.opt = 1
    self.My_Update(request_new_data = True)

  def OnRadio3(self, event = None):
    print "clicked on Radio3"

  def OnRadio4(self, event = None):
    print "clicked on Radio4"

  def OnSize(self, event = None):
    if( self.sizing_counter > 5 ):
      siz_data = self.GetSize()
      #print "New size of window =", siz_data
      optm_aspec_ratio = 4.21
      print "siz_data = ", siz_data[0], siz_data[1]
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
      self.My_Update(request_new_data = False)
      #print "resizing"
      print "aspec_ratio =", aspec_ratio
    else:
      self.sizing_counter += 1

  def My_Update(self, request_new_data = True):

    if( request_new_data == True ):
      self.dat, self.bkg, self.msk = self.tabl(opt = self.opt)
      self.I_max = self.tabl.Get_Max()
      self.box_lmt = self.tabl.Get_bbox()


      self.wx_Img_01, self.img_width, self.img_height = self.bmp(
                                      np_img_2d = self.dat
                                    , Intst_max = self.I_max
                                    , ofst = self.box_lmt)

      self.wx_Img_02, self.img_width, self.img_height = self.bmp(
                                      np_img_2d = self.bkg
                                    , Intst_max = self.I_max
                                    , ofst = self.box_lmt)

      self.wx_Img_03, self.img_width, self.img_height = self.bmp(
                                      np_img_2d = self.msk
                                    #, Intst_max = -1
                                    , Intst_max = self.I_max
                                    , ofst = self.box_lmt)

    self.My_Img_01 = from_wx_image_to_wx_bitmap(self.wx_Img_01
            , self.img_width, self.img_height, self.frame_scale)

    self.Image_01.SetBitmap(self.My_Img_01)

    self.My_Img_02 = from_wx_image_to_wx_bitmap(self.wx_Img_02
            , self.img_width, self.img_height, self.frame_scale)

    self.Image_02.SetBitmap(self.My_Img_02)

    self.My_Img_03 = from_wx_image_to_wx_bitmap(self.wx_Img_03
            , self.img_width, self.img_height, self.frame_scale)

    self.Image_03.SetBitmap(self.My_Img_03)

    self.Layout()
    self.Refresh()


  def On_Close_Window(self, event):
    self.Destroy()

class App(wx.App):
  def OnInit(self):
    self.frame = ReflectionFrame(None, -1, "DIALS Reflections Viewer"
    , wx.DefaultPosition,(550,200))

    self.SetTopWindow(self.frame)
    self.frame.Show(True)
    return True

  def table_in(self, loc_tabl):
    self.frame.tabl_to_frame(loc_tabl)


if __name__ == "__main__":

  import cPickle as pickle
  import sys
  from dials.model.data import ReflectionList # implicit import
  table = pickle.load(open(sys.argv[1]))
  pos_of_best_example_so_far = '''
  .... dials_regression/refinement_test_data/radiation_damaged_thaumatin/
  indexed.pickle
  '''
  #print "num of ref =", len(table)

  app = App(redirect=False)
  app.table_in(table)

  app.MainLoop()

