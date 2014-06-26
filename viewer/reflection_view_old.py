from __future__ import division
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

class MyFrame(wx.Frame):
  def __init__(self, *args, **kwargs):
    wx.Frame.__init__(self, *args, **kwargs)

    self.MaxImageSize = 300

    btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
    btn_prv_refl = wx.Button(self, -1, "Previous Reflection")
    btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext_refl)
    btn_prv_refl.Bind(wx.EVT_BUTTON, self.DisplayPrev_refl)

    btn_nxt_slice = wx.Button(self, -1, "Next slice ")
    btn_prv_slice = wx.Button(self, -1, "Previous slice")

    btn_chg_displ = wx.Button(self, -1, "change display")

    btn_nxt_slice.Bind(wx.EVT_BUTTON, self.DisplayNext_slice)
    btn_prv_slice.Bind(wx.EVT_BUTTON, self.DisplayPrev_slice)
    btn_chg_displ.Bind(wx.EVT_BUTTON, self.ChangeDisplay)
    # starting with an EmptyBitmap
    self.Image_01 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSize, self.MaxImageSize))

    # Using a Sizers to handle the layout

    v_box = wx.BoxSizer(wx.VERTICAL)
    h_box = wx.BoxSizer(wx.HORIZONTAL)
    u_box = wx.BoxSizer(wx.HORIZONTAL)

    u_box.Add(btn_prv_refl, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,5)

    v_box.Add(u_box)

    h_box.Add(self.Image_01
            , 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE
            , 7)
    r_box = wx.BoxSizer(wx.VERTICAL)
    r_box.Add(btn_nxt_slice, 0, wx.CENTER | wx.ALL,5)
    r_box.Add(btn_prv_slice, 0, wx.CENTER | wx.ALL,5)
    r_box.Add(btn_chg_displ, 0, wx.CENTER | wx.ALL,5)
    h_box.Add(r_box)
    v_box.Add(h_box)

    self.SetSizerAndFit(v_box)
    wx.EVT_CLOSE(self, self.OnCloseWindow)

  def tabl_to_frame(self, loc_tabl):
    self.tabl = table_s_navigator(loc_tabl)
    options_to_show = ''' bkg, dat, msk '''
    self.to_show = "dat"
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

  def ChangeDisplay(self, event = None):
    print "change display"
    options_to_show = ''' bkg, dat, msk '''
    if self.to_show == "dat":
      self.to_show = "bkg"
    elif self.to_show == "bkg":
      self.to_show = "msk"
    else:
      self.to_show = "dat"
    self.My_Update()

  def My_Update(self):
    bkg, dat, msk = self.tabl()
    if self.to_show == "dat":
      np_img = dat
    elif self.to_show == "bkg":
      np_img = bkg
    else:
      np_img = msk

    My_Img = GetBitmap_from_np_array(np_img)
    self.Image_01.SetBitmap(My_Img)
    self.Fit()
    self.Layout()
    self.Refresh()

  def OnCloseWindow(self, event):
    self.Destroy()

class App(wx.App):
  def OnInit(self):
    self.frame = MyFrame(None, -1, "DIALS Reflections Viewer"
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
