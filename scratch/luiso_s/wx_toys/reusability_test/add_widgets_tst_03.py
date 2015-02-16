#
#  DIALS viewer_frame test
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#  started from an example downloaded from: http://www.blog.pythonlibrary.org
#  then evolved to my needs
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

import wx
import wx.lib.scrolledpanel as scroll_pan

from bitmap_from_numpy \
     import GetBitmap_from_np_array, build_np_img

#class MyPanel(wx.Panel):
#    def __init__(self, parent):
#        wx.Panel.__init__(self, parent)


class multi_img_scrollable(scroll_pan.ScrolledPanel):

    def __init__(self, outer_panel):
        print "__init__ 01"
        super(multi_img_scrollable, self).__init__(outer_panel)


        self.number_of_img = 0
        self.frame = outer_panel

        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        controlSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.widgetSizer = wx.BoxSizer(wx.HORIZONTAL)

        self.addButton = wx.Button(self, label="Add")
        self.addButton.Bind(wx.EVT_BUTTON, self.onAddWidget)
        controlSizer.Add(self.addButton, 0, wx.CENTER|wx.ALL, 5)

        self.removeButton = wx.Button(self, label="Remove")
        self.removeButton.Bind(wx.EVT_BUTTON, self.onRemoveWidget)
        controlSizer.Add(self.removeButton, 0, wx.CENTER|wx.ALL, 5)

        self.mainSizer.Add(controlSizer, 0, wx.CENTER)
        self.mainSizer.Add(self.widgetSizer, 0, wx.CENTER|wx.ALL, 10)

        self.SetSizer(self.mainSizer)


        self.SetupScrolling()
        print "__init__ 02"


    def onAddWidget(self, event):
        self.number_of_img += 1
        label = "Button %s" %  self.number_of_img
        name = "button%s" % self.number_of_img


        data2d = build_np_img(width=5, height=8)
        bitmap = GetBitmap_from_np_array(data2d)
        bitmap_tmp = wx.StaticBitmap(self, bitmap=bitmap)
        self.widgetSizer.Add(bitmap_tmp, 0, wx.ALL, 5)

        self.frame.fSizer.Layout()
        #self.frame.Fit()
        #self.Refresh()


    def onRemoveWidget(self, event):
        if self.widgetSizer.GetChildren():
            self.widgetSizer.Hide(self.number_of_img-1)
            self.widgetSizer.Remove(self.number_of_img-1)
            self.number_of_img -= 1
            self.frame.fSizer.Layout()
            #self.frame.Fit()
            #self.Refresh()
            print "number_of_img =", self.number_of_img

class MyFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, parent=None, title="Add / Remove Buttons")
        self.fSizer = wx.BoxSizer(wx.VERTICAL)
        panel = multi_img_scrollable(self)

        self.fSizer.Add(panel, 1, wx.EXPAND)
        self.fSizer.SetMinSize((400, 200))
        self.SetSizer(self.fSizer)
        self.Fit()
        self.Show()

#----------------------------------------------------------------------
if __name__ == "__main__":
    app = wx.App(False)
    frame = MyFrame()
    app.MainLoop()
