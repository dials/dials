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
from dials.array_family import flex

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
    btn_nxt_slice.Bind(wx.EVT_BUTTON, self.DisplayNext_slice)
    btn_prv_slice.Bind(wx.EVT_BUTTON, self.DisplayPrev_slice)

    # starting with an EmptyBitmap
    self.Image = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSize, self.MaxImageSize))

    self.DisplayNext_refl()

    # Using a Sizers to handle the layout

    v_box = wx.BoxSizer(wx.VERTICAL)
    h_box = wx.BoxSizer(wx.HORIZONTAL)
    u_box = wx.BoxSizer(wx.HORIZONTAL)

    u_box.Add(btn_prv_refl, 0, wx.CENTER | wx.ALL,5)
    u_box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,5)

    v_box.Add(u_box)

    h_box.Add(self.Image
            , 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE
            , 7)
    r_box = wx.BoxSizer(wx.VERTICAL)
    r_box.Add(btn_nxt_slice, 0, wx.CENTER | wx.ALL,5)

    r_box.Add(btn_prv_slice, 0, wx.CENTER | wx.ALL,5)
    h_box.Add(r_box)
    v_box.Add(h_box)

    self.SetSizerAndFit(v_box)
    wx.EVT_CLOSE(self, self.OnCloseWindow)

    #self.tmp_img = 0

  def DisplayNext_refl(self, event = None):
    np_img = build_np_img(width = 20, height = 30)
    self.My_Update(np_img)

  def DisplayPrev_refl(self, event = None):
    np_img = build_np_img(width = 20, height = 10)
    self.My_Update(np_img)

  def DisplayNext_slice(self, event = None):
    #np_img = build_np_img(width = 25, height = 15)
    #self.My_Update(np_img)
    self.My_Update(self.tmp_img)

  def DisplayPrev_slice(self, event = None):
    np_img = build_np_img(width = 10, height = 15)
    self.My_Update(np_img)

  def My_Update(self, np_img):
    My_Img = GetBitmap_from_np_array(np_img)
    self.Image.SetBitmap(My_Img)
    self.Fit()
    self.Layout()
    self.Refresh()

  def OnCloseWindow(self, event):
    self.Destroy()

class App(wx.App):
  def OnInit(self):
    frame = MyFrame(None, -1, "DIALS Reflections Viewer"
    , wx.DefaultPosition,(550,200))
    self.SetTopWindow(frame)
    frame.Show(True)
    return True

class table_s_navigator(object):

  def __init__(self, table):

    self.table = table
    self.num_ref = len(table)
    if self.num_ref >= 1:
      row = table[52]
    else:
      print "ERROR 0 reflections"

    self.data_flex = row['shoebox'].data
    self.background_flex = row['shoebox'].background
    self.mask_flex = row['shoebox'].mask
    self.depth = self.data_flex.all()[0]
    print "depth of refl =", self.depth
    if self.depth >= 1:
      self.z = 0
    else:
      print "ERROR 0 depth"

  def next_slice(self):
    if self.z < self.depth:
      self.z += 1
      self.__call__()
    else:
      print "maximum depth reached"
  def Previous_slice(self):
    pass
  def next_Reflection(self):
    pass
  def Previous_Reflection(self):
    pass

  def __call__(self):

    mask2d = self.mask_flex[self.z:self.z + 1, :, :]
    mask2d.reshape(flex.grid(self.mask_flex.all()[1:]))
    mask2d_np = mask2d.as_numpy_array()

    data2d = self.data_flex[self.z:self.z + 1, :, :]
    data2d.reshape(flex.grid(self.data_flex.all()[1:]))
    data2d_np= data2d.as_numpy_array()

    background2d = self.background_flex[self.z:self.z + 1, :, :]
    background2d.reshape(flex.grid(self.background_flex.all()[1:]))
    background2d_np = background2d.as_numpy_array()

    self.img_background = background2d_np
    self.img_data = data2d_np
    self.img_mask = mask2d_np
    return self.img_background, self.img_data, self.img_mask

  def background(self):
    print "from background(self)"
    return self.img_background

  def data(self):
    print "from data(self)"
    return self.img_data

  def mask(self):
    print "from mask(self)"
    return self.img_mask

old_mod = '''
def paint_refl(table):
  app = App(redirect=False)

  tbl = table_s_navigator(table)
  bkg, dat, msk = tbl()
  MyFrame.tmp_img = msk
  #MyFrame.tmp_img = tbl.data()
  app.MainLoop()

  return
#'''
if __name__ == "__main__":

  import cPickle as pickle
  import sys
  from dials.model.data import ReflectionList # implicit import
  table = pickle.load(open(sys.argv[1]))

  print "num of ref =", len(table)

  #paint_refl(table)
  app = App(redirect=False)
  
  tbl = table_s_navigator(table)
  bkg, dat, msk = tbl()
  MyFrame.tmp_img = dat  

  app.MainLoop()
#'''