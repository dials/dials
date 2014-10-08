import os
import wx
import numpy as np
import wx.lib.scrolledpanel as scroll_pan
from bitmap_from_numpy import GetBitmap_from_np_array

class multi_img(scroll_pan.ScrolledPanel):
    def __init__(self, outer_frame):
        scroll_pan.ScrolledPanel.__init__(self, outer_frame)

        print "Hi from scroll_pan"

        data2d = np.arange( 7 * 6, dtype = 'uintc').reshape( 7, 6 )
        bmp = GetBitmap_from_np_array(data2d)
        local_bmp = wx.StaticBitmap(self, bitmap=bmp)

        data2d_ag = np.arange( 6 * 7, dtype = 'uintc').reshape( 6, 7 )
        bmp_ag = GetBitmap_from_np_array(data2d_ag)
        local_bmp_ag = wx.StaticBitmap(self, bitmap=bmp_ag)

        self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.top_sizer.Add(local_bmp, 1, wx.RIGHT)
        self.top_sizer.Add(local_bmp_ag, 1, wx.LEFT)
        self.SetSizer(self.top_sizer)
        self.SetupScrolling()

class iner_panel(wx.Panel):
    def __init__(self, outer_frame):
        wx.Panel.__init__(self, outer_frame)

        print "Hi"

        data2d = np.arange( 7 * 6, dtype = 'uintc').reshape( 7, 6 )
        bmp = GetBitmap_from_np_array(data2d)
        local_bmp = wx.StaticBitmap(self, bitmap=bmp)

        data2d_ag = np.arange( 6 * 7, dtype = 'uintc').reshape( 6, 7 )
        bmp_ag = GetBitmap_from_np_array(data2d_ag)
        local_bmp_ag = wx.StaticBitmap(self, bitmap=bmp_ag)

        self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.top_sizer.Add(local_bmp, 1, wx.RIGHT)
        self.top_sizer.Add(local_bmp_ag, 1, wx.LEFT)
        self.SetSizer(self.top_sizer)

class MyFrame(wx.Frame):
    """ We simply derive a new class of Frame. """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=wx.DefaultSize)

        # Attributes
        self.panel_01 = iner_panel(self)
        self.panel_02 = multi_img(self)
        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.panel_01, 1, wx.EXPAND)
        sizer.Add(self.panel_02, 1, wx.EXPAND)
        self.SetSizer(sizer)

        self.Show(True)
app = wx.App(False)
frame = MyFrame(None, 'Test reuse')
app.MainLoop()
