import os
import wx
import numpy as np

from bitmap_from_numpy import GetBitmap_from_np_array

class iner_panel(wx.Panel):
    def __init__(self, outer_frame):
        wx.Panel.__init__(self, outer_frame)

        print "Hi"

        data2d = np.arange( 5 * 5, dtype = 'uintc').reshape( 5, 5 )
        bmp = GetBitmap_from_np_array(data2d)
        local_bmp = wx.StaticBitmap(self, bitmap=bmp)

        #to_test = '''
        self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.top_sizer.Add(local_bmp, 1, wx.EXPAND)
        self.SetSizer(self.top_sizer)
        #'''

class MyFrame(wx.Frame):
    """ We simply derive a new class of Frame. """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=wx.DefaultSize)

        # Attributes
        self.panel = iner_panel(self)

        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.panel, 1, wx.EXPAND)
        self.SetSizer(sizer)

        self.Show(True)
app = wx.App(False)
frame = MyFrame(None, 'Test reuse')
app.MainLoop()
