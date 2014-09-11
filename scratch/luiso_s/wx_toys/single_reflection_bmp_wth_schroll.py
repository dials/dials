#
#  DIALS viewer_frame
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."
#
import wx
import wx.lib.scrolledpanel as scroll_pan

from dials.viewer.img_utilities import GetBitmap_from_np_array

from dials.scratch.luiso_s.wx_toys.bitmap_from_numpy_w_matplotlib_well_done \
     import build_np_img
class ImageListCtrl(scroll_pan.ScrolledPanel):
    """Simple control to display a list of images"""
    def __init__(self, parent, bitmaps=list(),
                 style=wx.TAB_TRAVERSAL|wx.BORDER_SUNKEN):
        super(ImageListCtrl, self).__init__(parent,
                                            style=style)

        # Attributes
        self.images = list()
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        # Setup
        for bmp in bitmaps:
            self.AppendBitmap(bmp)
        self.SetSizer(self.sizer)

    def AppendBitmap(self, bmp):
        """Add another bitmap to the control"""
        self.images.append(bmp)
        sbmp = wx.StaticBitmap(self, bitmap=bmp)
        self.sizer.Add(sbmp, 0, wx.EXPAND|wx.TOP, 5)
        self.SetupScrolling()


class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyFrame(None, title="ScrolledPanel", size=(300,200))
        self.SetTopWindow(self.frame)
        self.frame.Show()

        return True

class MyFrame(wx.Frame):
    def __init__(self, parent, *args, **kwargs):
        wx.Frame.__init__(self, parent, *args, **kwargs)

        # Attributes
        self.panel = MyPanel(self)

        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.panel, 1, wx.EXPAND)
        self.SetSizer(sizer)

class MyPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

        # Attributes
        self.il = ImageListCtrl(self)

        # Setup
        for times in range(5):
            data2d = build_np_img(width=5, height=3)
            my_bitmap = GetBitmap_from_np_array(data2d)
            self.il.AppendBitmap(my_bitmap)

        # Layout
        self.mid_sizer = wx.BoxSizer(wx.VERTICAL)
        self.mid_sizer.Add(self.il, 1, wx.EXPAND)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(wx.StaticText(self, label="Image List:"), 0)
        self.sizer.Add(self.mid_sizer, 1, wx.EXPAND)

        self.SetSizer(self.sizer)

if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
