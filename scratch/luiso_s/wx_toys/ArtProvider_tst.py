#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."
#
#  learning from an example given by Andrea Gavana
#  This sample shows how to retrieve default platform's
#  bitmaps using wx.ArtProvider

import wx

class BitmapFrame(wx.Frame):

    def __init__(self):

        wx.Frame.__init__(self, None, -1, title='ArtProvider example')

        panel = wx.Panel(self)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        bitmap_sizer = wx.BoxSizer(wx.HORIZONTAL)

        bitmap_sizer.Add((0, 0), 1, wx.EXPAND)

        # Show a few bitmaps retrieved via wx.ArtProvider
        for kind in [wx.ART_GO_UP, wx.ART_GO_DOWN, wx.ART_GO_BACK, wx.ART_GO_FORWARD]:
            bmp = wx.ArtProvider.GetBitmap(kind, wx.ART_OTHER, (32, 32))
            static_bitmap = wx.StaticBitmap(panel, -1, bmp)
            bitmap_sizer.Add(static_bitmap, 0, wx.ALL, 5)

        # Layout everything in a nice sizer
        bitmap_sizer.Add((0, 0), 1, wx.EXPAND)
        main_sizer.Add((0, 0), 1, wx.EXPAND)
        main_sizer.Add(bitmap_sizer, 0, wx.EXPAND)
        main_sizer.Add((0, 0), 1, wx.EXPAND)

        panel.SetSizer(main_sizer)
        main_sizer.SetSizeHints(panel)
        main_sizer.Layout()


if __name__ == '__main__':
    app = wx.App(0)
    frame = BitmapFrame()
    frame.Show()
    app.MainLoop()