#
#  DIALS viewer_frame
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#                   and
#          "wxPython aplication Development CookBook" authors
#
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

#
import wx
import wx.lib.scrolledpanel as scrolledpanel

from dials.scratch.luiso_s.wx_toys.bitmap_from_numpy_w_matplotlib_well_done \
     import GetBitmap_from_np_array, build_np_img
class ImageListCtrl(scrolledpanel.ScrolledPanel):
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
        for art in (wx.ART_ERROR, wx.ART_WARNING, wx.ART_INFORMATION,
                    wx.ART_COPY, wx.ART_PASTE, wx.ART_CUT, wx.ART_CDROM,
                    wx.ART_HARDDISK, wx.ART_FOLDER, wx.ART_FLOPPY):
            bmp = wx.ArtProvider.GetBitmap(art)
            self.il.AppendBitmap(bmp)

        # Layout

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(wx.StaticText(self, label="Image List:"), 0)
        self.mid_sizer = wx.BoxSizer(wx.VERTICAL)
        self.mid_sizer.Add(self.il, 1, wx.EXPAND)
        self.sizer.Add(self.mid_sizer, 1, wx.EXPAND)

        self.addButton = wx.Button(self, label="new data")
        self.addButton.Bind(wx.EVT_BUTTON, self.OnNew_data)
        self.sizer.Add(self.addButton, 0, wx.CENTER|wx.ALL, 5)

        self.SetSizer(self.sizer)

    def OnNew_data(self, event):
        print "from OnNew_data"
        self.il = ImageListCtrl(self)

        print "here 01"
        self.mid_sizer.Hide(0)
        self.mid_sizer.Remove(0)
        print "here 02"

        '''
        for art in (wx.ART_ERROR, wx.ART_WARNING, wx.ART_INFORMATION,
                    wx.ART_COPY, wx.ART_PASTE, wx.ART_CUT, wx.ART_CDROM,
                    wx.ART_HARDDISK, wx.ART_FOLDER, wx.ART_FLOPPY):
            bmp = wx.ArtProvider.GetBitmap(art)
            self.il.AppendBitmap(bmp)
        '''

        for times in range(5):
            data2d = build_np_img(width=5, height=3)
            my_bitmap = GetBitmap_from_np_array(data2d)
            #bitmap_tmp = wx.StaticBitmap(self, bitmap=my_bitmap)
            self.il.AppendBitmap(my_bitmap)

        self.mid_sizer = wx.BoxSizer(wx.VERTICAL)
        self.mid_sizer.Add(self.il, 1, wx.EXPAND)

        print "here 03"
        self.sizer.Add(self.il, 1, wx.EXPAND)
        print "here 04"

        self.Refresh()
        self.Update()
        self.Layout()
        self.Fit()

        wx.Yield()


if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
