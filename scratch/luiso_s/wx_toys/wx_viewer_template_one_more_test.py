#
#  DIALS viewer
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

import wx, os
from dials.viewer.viewer_utilities \
     import GetBitmap_from_np_array, build_np_img, from_wx_image_to_wx_bitmap

import numpy

class ReflectionFrame(wx.Frame):
  def __init__(self, *args, **kwargs):
    wx.Frame.__init__(self, *args, **kwargs)


    self.MaxImageSizeX = 320
    self.MaxImageSizeY = 240

    btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
    btn_prv_refl = wx.Button(self, -1, "Previous Reflection")
    btn_nxt_slice = wx.Button(self, -1, "Next slice ")
    btn_prv_slice = wx.Button(self, -1, "Previous slice")
    btn_tst = wx.Button(self, -1, "tst btn")
    btn_tst1 = wx.Button(self, -1, "tst btn1")

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

    r_box.Add(btn_tst, 0, wx.CENTER | wx.ALL,5)
    r_box.Add(btn_tst1, 0, wx.CENTER | wx.ALL,5)

    h_box.Add(r_box)
    v_box.Add(h_box)

    self.frame_scale = 0.5

    self.sizing_counter = 0
    self.SetSizerAndFit(v_box)

    btn_tst.Bind(wx.EVT_BUTTON, self.B_tst)
    btn_tst1.Bind(wx.EVT_BUTTON, self.B_tst1)
    self.Bind(wx.EVT_SIZE, self.OnSize)

    wx.EVT_CLOSE(self, self.On_Close_Window)

  def tabl_to_frame(self):
    self.My_Update()
  def B_tst(self, event = None):
    self.frame_scale = self.frame_scale * 1.1
    self.My_Update()
    print "self.GetSize() =", self.GetSize()
  def B_tst1(self, event = None):
    self.frame_scale = self.frame_scale * 0.9
    self.My_Update()
    print "self.GetSize() =", self.GetSize()

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
      self.My_Update(request_new_size = False)
      print "resizing"
    else:
      self.sizing_counter += 1

  def My_Update(self, request_new_size = True):
    if( request_new_size == True ):
      np_tmp = numpy.zeros( (20, 20), 'float')
      self.bkg = np_tmp
      self.dat = np_tmp
      self.msk = np_tmp
      self.I_max = 100
      print "re - fitting"

      self.wx_Img_dat, self.img_width, self.img_height = GetBitmap_from_np_array(
                                      np_img_2d = self.dat
                                    , Intst_max = self.I_max
                                    , img_scale = self.frame_scale)

      self.wx_Img_bkg, self.img_width, self.img_height = GetBitmap_from_np_array(
                                      np_img_2d = self.bkg
                                    , Intst_max = self.I_max
                                    , img_scale = self.frame_scale)

      self.wx_Img_msk, self.img_width, self.img_height = GetBitmap_from_np_array(
                                      np_img_2d = self.msk
                                    , Intst_max = self.I_max
                                    , img_scale = self.frame_scale)

      self.Fit()
    else:
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

  def table_in(self):
    self.frame.tabl_to_frame()


if __name__ == "__main__":

  app = App(redirect=False)
  app.table_in()

  app.MainLoop()

