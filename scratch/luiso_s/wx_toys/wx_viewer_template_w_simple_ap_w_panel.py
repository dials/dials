#from dials.viewer.viewer_utilities import GetBitmap_from_np_array, build_np_img

from bitmap_from_numpy import GetBitmap_from_np_array

import wx

class RandomPanel(wx.Panel):
  def __init__(self, parent, style=wx.FULL_REPAINT_ON_RESIZE):
    wx.Panel.__init__(self, parent)

    self.MaxImageSizeX = 320
    self.MaxImageSizeY = 240

    btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
    btn_prv_refl = wx.Button(self, -1, "Previous Reflection")
    #btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext_refl)
    #btn_prv_refl.Bind(wx.EVT_BUTTON, self.DisplayPrev_refl)

    btn_nxt_slice = wx.Button(self, -1, "Next slice ")
    btn_prv_slice = wx.Button(self, -1, "Previous slice")
    btn_tst = wx.Button(self, -1, "tst btn")
    self.Image_01 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_02 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))
    self.Image_03 = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                 self.MaxImageSizeX, self.MaxImageSizeY))

    #btn_nxt_slice.Bind(wx.EVT_BUTTON, self.DisplayNext_slice)
    #btn_prv_slice.Bind(wx.EVT_BUTTON, self.DisplayPrev_slice)

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

    #self.is_sizing_flag = False
    self.frame_scale = 0.5

    self.SetSizerAndFit(v_box)
    btn_tst.Bind(wx.EVT_BUTTON, self.tst_tmp)
    self.Bind(wx.EVT_SIZE, self.OnSize)
    self.sizing_counter = 1
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

  def tst_tmp(self, event = None):
    self.My_Update()

  def My_Update(self, request_new_size = True):
    image_np = build_np_img(20, 20)
    wxBitmap = GetBitmap_from_np_array(image_np, 100, self.frame_scale)
    self.Image_01 = wxBitmap
    self.Image_02 = wxBitmap
    self.Image_03 = wxBitmap

    if( request_new_size == True ):
      self.Fit()

    self.Layout()
    self.Refresh()



if __name__ == "__main__":
  app = wx.PySimpleApp()
  frame = wx.Frame(None, title='panel with reflection viewer')
  panel = RandomPanel(frame)
  #panel.draw()
  frame.Show()
  app.MainLoop()
