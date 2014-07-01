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
     import GetBitmap_from_np_array, from_wx_image_to_wx_bitmap
#from dials.viewer.viewer_utilities import build_np_img
from dials.viewer.reflection_data_navigator import table_s_navigator


class ReflectionFrame(wx.Frame):
  def __init__(self, *args, **kwargs):
    wx.Frame.__init__(self, *args, **kwargs)


    self.MaxImageSizeX = 220
    self.MaxImageSizeY = 140

    btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
    btn_prv_refl = wx.Button(self, -1, "Previous Reflection")
    btn_nxt_slice = wx.Button(self, -1, "Next slice ")
    btn_prv_slice = wx.Button(self, -1, "Previous slice")

    test_code = '''
    btn_tst = wx.Button(self, -1, "tst btn")
    btn_tst1 = wx.Button(self, -1, "tst btn1")
    '''

    self.Image_01 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_02 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_03 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))

    v_box = wx.BoxSizer(wx.VERTICAL)
    h_box = wx.BoxSizer(wx.HORIZONTAL)
    u_box = wx.BoxSizer(wx.HORIZONTAL)

    u_box.Add(btn_prv_refl, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,5)

    v_box.Add(u_box)

    h_box.Add(self.Image_01, 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 7)

    h_box.Add(self.Image_02, 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 7)

    h_box.Add(self.Image_03, 0
            , wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 7)

    r_box = wx.BoxSizer(wx.VERTICAL)
    r_box.Add(btn_nxt_slice, 0, wx.CENTER | wx.ALL,5)
    r_box.Add(btn_prv_slice, 0, wx.CENTER | wx.ALL,5)
    r_box.AddSpacer(50)

    test_code = '''
    r_box.Add(btn_tst, 0, wx.CENTER | wx.ALL,5)
    r_box.Add(btn_tst1, 0, wx.CENTER | wx.ALL,5)
    '''

    h_box.Add(r_box)
    v_box.Add(h_box)

    self.frame_scale = 0.3

    self.sizing_counter = 0
    self.SetSizerAndFit(v_box)

    btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext_refl)
    btn_prv_refl.Bind(wx.EVT_BUTTON, self.DisplayPrev_refl)
    btn_nxt_slice.Bind(wx.EVT_BUTTON, self.DisplayNext_slice)
    btn_prv_slice.Bind(wx.EVT_BUTTON, self.DisplayPrev_slice)
    test_code = '''
    btn_tst.Bind(wx.EVT_BUTTON, self.B_tst)
    btn_tst1.Bind(wx.EVT_BUTTON, self.B_tst1)
    '''
    self.Bind(wx.EVT_SIZE, self.OnSize)

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
  test_code = '''
  def B_tst(self, event = None):
    self.frame_scale = self.frame_scale * 1.1
    self.My_Update()
    print "self.GetSize() =", self.GetSize()
  def B_tst1(self, event = None):
    self.frame_scale = self.frame_scale * 0.9
    self.My_Update()
    print "self.GetSize() =", self.GetSize()
  '''

  def OnSize(self, event = None):
    if( self.sizing_counter > 5 ):
      siz_data = self.GetSize()
      print "New size of window =", siz_data
      optm_aspec_ratio = 3.81
      print "siz_data = ", siz_data[0], siz_data[1]
      print "aspect ratio = ", float(siz_data[0])/ float(siz_data[1])
      aspec_ratio = float(siz_data[0])/ float(siz_data[1])

      if(aspec_ratio > optm_aspec_ratio):
        print "use float(siz_data[1] (Height) to calculate new size"
        self.frame_scale = float(siz_data[1]) * 0.5 / 291.0
        #(1100, 291)
      else:
        print "use float(siz_data[0] (with) to calculate new size"
        self.frame_scale = float(siz_data[0]) * 0.5 / 1100.0
        #(1100, 291)
      self.My_Update(request_new_data = False)
      print "resizing"
    else:
      self.sizing_counter += 1

  def My_Update(self, request_new_data = True):

    if( request_new_data == True ):
      self.bkg, self.dat, self.msk = self.tabl()
      self.I_max = self.tabl.Get_Max()
      self.box_lmt = self.tabl.Get_bbox()

      self.wx_Img_dat, self.img_width, self.img_height = GetBitmap_from_np_array(
                                      np_img_2d = self.dat
                                    , Intst_max = self.I_max
                                    , img_scale = self.frame_scale
                                    , ofst = self.box_lmt)

      self.wx_Img_bkg, self.img_width, self.img_height = GetBitmap_from_np_array(
                                      np_img_2d = self.bkg
                                    , Intst_max = self.I_max
                                    , img_scale = self.frame_scale
                                    , ofst = self.box_lmt)

      self.wx_Img_msk, self.img_width, self.img_height = GetBitmap_from_np_array(
                                      np_img_2d = self.msk
                                    , Intst_max = -1
                                    , img_scale = self.frame_scale
                                    , ofst = self.box_lmt)

    self.My_Img_01 = from_wx_image_to_wx_bitmap(self.wx_Img_dat
            , self.img_width, self.img_height, self.frame_scale)

    self.Image_01.SetBitmap(self.My_Img_01)

    self.My_Img_02 = from_wx_image_to_wx_bitmap(self.wx_Img_bkg
            , self.img_width, self.img_height, self.frame_scale)

    self.Image_02.SetBitmap(self.My_Img_02)

    self.My_Img_03 = from_wx_image_to_wx_bitmap(self.wx_Img_msk
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
  print "num of ref =", len(table)

  app = App(redirect=False)
  app.table_in(table)

  app.MainLoop()

