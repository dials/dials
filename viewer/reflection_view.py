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
from dials.viewer.viewer_utilities import GetBitmap_from_np_array, build_np_img
from dials.viewer.reflection_data_navigator import table_s_navigator

class ReflectionFrame(wx.Frame):
  def __init__(self, *args, **kwargs):
    wx.Frame.__init__(self, *args, **kwargs)

    self.MaxImageSizeX = 320
    self.MaxImageSizeY = 240

    btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
    btn_prv_refl = wx.Button(self, -1, "Previous Reflection")
    btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext_refl)
    btn_prv_refl.Bind(wx.EVT_BUTTON, self.DisplayPrev_refl)

    btn_nxt_slice = wx.Button(self, -1, "Next slice ")
    btn_prv_slice = wx.Button(self, -1, "Previous slice")
    btn_tst = wx.Button(self, -1, "tst btn")

    btn_nxt_slice.Bind(wx.EVT_BUTTON, self.DisplayNext_slice)
    btn_prv_slice.Bind(wx.EVT_BUTTON, self.DisplayPrev_slice)
    btn_tst.Bind(wx.EVT_BUTTON, self.B_tst)

    self.Image_01 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_02 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_03 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Bind(wx.EVT_SIZE, self.OnSize)
    # Using a Sizers to handle the layout

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

    h_box.Add(r_box)
    v_box.Add(h_box)

    self.frame_scale = 0.5

    self.SetSizerAndFit(v_box)
    wx.EVT_CLOSE(self, self.OnCloseWindow)

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
  def B_tst(self, event = None):
    print "Here tst"
    self.frame_scale = self.frame_scale * 1.2
    self.My_Update()
    print "self.GetSize() =", self.GetSize()

  def OnSize(self, event = None):
    siz_data = self.GetSize()
    print "New size of window =", siz_data
    #self.frame_scale = float(siz_data[0] * siz_data[1]) * 0.5 / 320100.0
    #self.My_Update()

  def My_Update(self):
    #self.Unbind(wx.EVT_SIZE)
    bkg, dat, msk = self.tabl()

    I_max = self.tabl.Get_Max()

    My_Img = GetBitmap_from_np_array(np_img_2d = dat, Intst_max = I_max
                                     , img_scale = self.frame_scale)
    self.Image_01.SetBitmap(My_Img)
    My_Img = GetBitmap_from_np_array(np_img_2d = bkg, Intst_max = I_max
                                     , img_scale = self.frame_scale)
    self.Image_02.SetBitmap(My_Img)
    My_Img = GetBitmap_from_np_array(np_img_2d = msk, Intst_max = -1
                                     , img_scale = self.frame_scale)
    self.Image_03.SetBitmap(My_Img)
    self.Fit()
    self.Layout()
    self.Refresh()

  def OnCloseWindow(self, event):
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

  print "num of ref =", len(table)

  app = App(redirect=False)
  app.table_in(table)

  app.MainLoop()
#'''
