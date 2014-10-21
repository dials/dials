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

        my_sizer = wx.BoxSizer(wx.HORIZONTAL)
        my_sizer.Add(local_bmp, 10, wx.RIGHT)
        my_sizer.Add(local_bmp_ag, 10, wx.LEFT)
        self.SetSizer(my_sizer)
        self.SetupScrolling()

class butoms_panel(wx.Panel):
    def __init__(self, outer_frame):
        wx.Panel.__init__(self, outer_frame)

        Hide_I_Button = wx.Button(self, label="Hide I")
        Hide_I_Button.Bind(wx.EVT_BUTTON, self.OnShwIBut)

        Show_I_Button = wx.Button(self, label="Show I")
        Show_I_Button.Bind(wx.EVT_BUTTON, self.OnHidIBut)
        sb = wx.StaticBox(self, -1, size=(180, 150))
        my_sizer = wx.StaticBoxSizer(sb, wx.VERTICAL)
        my_sizer.Add(Show_I_Button, 0, wx.LEFT | wx.ALL,8)
        my_sizer.Add(Hide_I_Button, 0, wx.LEFT | wx.ALL,8)
        my_sizer.SetMinSize((120, 250))
        #my_sizer.SetMaxSize((140, 270))

        self.SetSizer(my_sizer)

    def OnShwIBut(self, event):
        print "OnShwIBut"

    def OnHidIBut(self, event):
        print "OnHidIBut"

class My_new_Frame(wx.Frame):
    """ We simply derive a new class of Frame. """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=wx.DefaultSize)

        # Attributes
        self.panel_01 = butoms_panel(self)
        self.panel_02 = multi_img(self)
        # Layout
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.panel_01, 1, wx.EXPAND)
        sizer.Add(self.panel_02, 1, wx.EXPAND)
        self.SetSizerAndFit(sizer)

        self.Show(True)
app = wx.App(False)
frame = My_new_Frame(None, 'Test reuse')
app.MainLoop()
