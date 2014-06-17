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

class TestFrame(wx.Frame):
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

        # Using a Sizer to handle the layout

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

    def DisplayNext_refl(self, event = None):
        np_img = build_np_img(width = 20, height = 30)
        self.My_Update(np_img)

    def DisplayPrev_refl(self, event = None):
        np_img = build_np_img(width = 20, height = 10)
        self.My_Update(np_img)

    def DisplayNext_slice(self, event = None):
        np_img = build_np_img(width = 25, height = 15)
        self.My_Update(np_img)

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

        frame = TestFrame(None, -1, "wxBitmap Test", wx.DefaultPosition
        ,(550,200))
        self.SetTopWindow(frame)
        frame.Show(True)
        return True

if __name__ == "__main__":
    app = App(0)
    app.MainLoop()