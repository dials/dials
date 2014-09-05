#
#  DIALS viewer_frame
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

from dials.scratch.luiso_s.wx_toys.bitmap_from_numpy_w_matplotlib_well_done \
     import GetBitmap_from_np_array, build_np_img

class MyPanel(wx.Panel):

    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.number_of_img = 0
        self.number_of_floors = 1
        self.frame = parent

        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        controlSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.widgetSizer = [wx.BoxSizer(wx.HORIZONTAL)]

        self.addButton = wx.Button(self, label="Add")
        self.addButton.Bind(wx.EVT_BUTTON, self.onAddWidget)
        controlSizer.Add(self.addButton, 0, wx.CENTER|wx.ALL, 5)

        self.removeButton = wx.Button(self, label="Remove")
        self.removeButton.Bind(wx.EVT_BUTTON, self.onRemoveWidget)
        controlSizer.Add(self.removeButton, 0, wx.CENTER|wx.ALL, 5)

        self.V_addButton = wx.Button(self, label="Vertical add")
        self.V_addButton.Bind(wx.EVT_BUTTON, self.on_V_add)
        controlSizer.Add(self.V_addButton, 0, wx.CENTER|wx.ALL, 5)

        self.V_rmButton = wx.Button(self, label="Vertical remove")
        self.V_rmButton.Bind(wx.EVT_BUTTON, self.on_V_rm)
        controlSizer.Add(self.V_rmButton, 0, wx.CENTER|wx.ALL, 5)



        #test button

        self.test_mem = wx.Button(self, label="test")
        self.test_mem.Bind(wx.EVT_BUTTON, self.on_test_mem)
        controlSizer.Add(self.test_mem, 0, wx.CENTER|wx.ALL, 5)


        #

        self.mainSizer.Add(controlSizer, 0, wx.CENTER)
        self.mainSizer.Add(self.widgetSizer[0], 0, wx.CENTER|wx.ALL, 10)

        self.SetSizer(self.mainSizer)

    def onAddWidget(self, event):
        self.number_of_img += 1
        label = "Button %s" %  self.number_of_img
        name = "button%s" % self.number_of_img

        for floor in range(self.number_of_floors):
            data2d = build_np_img(width=5, height=3)
            bitmap = GetBitmap_from_np_array(data2d)
            bitmap_tmp = wx.StaticBitmap(self, bitmap=bitmap)
            self.widgetSizer[floor].Add(bitmap_tmp, 0, wx.ALL, 5)

        self.frame.fSizer.Layout()
        self.frame.Fit()

    def onRemoveWidget(self, event):
        if self.widgetSizer[0].GetChildren():
            for floor in range(self.number_of_floors):
                self.widgetSizer[floor].Hide(self.number_of_img-1)
                self.widgetSizer[floor].Remove(self.number_of_img-1)

            self.number_of_img -= 1
            self.frame.fSizer.Layout()
            self.frame.Fit()
            print "number_of_img =", self.number_of_img

    def on_V_add(self, event):
        print "from on_V_add"
        self.widgetSizer.append(wx.BoxSizer(wx.HORIZONTAL))
        self.mainSizer.Add(self.widgetSizer[self.number_of_floors], 0, wx.CENTER|wx.ALL, 10)

        for n_siz in range(self.number_of_img):
            data2d = build_np_img(width=5, height=3)
            bitmap = GetBitmap_from_np_array(data2d)
            bitmap_tmp = wx.StaticBitmap(self, bitmap=bitmap)
            self.widgetSizer[self.number_of_floors].Add(bitmap_tmp, 0, wx.ALL, 5)

        self.number_of_floors += 1
        self.frame.fSizer.Layout()
        self.frame.Fit()
    def on_V_rm(self, event):

        floor = self.number_of_floors - 1

        '''
        for floor in range(self.number_of_floors):
                self.widgetSizer[floor].Hide(self.number_of_img-1)
                self.widgetSizer[floor].Remove(self.number_of_img-1)
        '''

        for n_siz in range(self.number_of_img -1, -1, -1):

            print "n_siz =", n_siz
            self.widgetSizer[floor].Hide(n_siz)
            self.widgetSizer[floor].Remove(n_siz)

        self.mainSizer.Hide(self.widgetSizer[floor])
        self.mainSizer.Remove(self.widgetSizer[floor])
        del self.widgetSizer[floor]

        self.number_of_floors = floor

        print "from V_rmButton"

        self.frame.fSizer.Layout()
        self.frame.Fit()


    def on_test_mem(self, event):
        for times_times in range(50):
            for times in range(4):
                self.onAddWidget(event)
            self.tst_update()

            for times in range(3):
                self.on_V_add(event)
            self.tst_update()

            for times in range(4):
                self.onRemoveWidget(event)
            self.tst_update()

            for times in range(3):
                self.on_V_rm(event)
            self.tst_update()

        print "test Done"

    def tst_update(self):
        import time
        # Wait for 5 seconds
        #time.sleep(1)
        self.frame.fSizer.Layout()
        #self.frame.Fit()
        self.frame.Refresh()
        self.frame.Update()
        wx.Yield()
        time.sleep(1)
        self.Layout()
        #self.Fit()
        self.Refresh()
        self.Update()
        wx.Yield()
        time.sleep(1)

        #self.Layout()

        remember_to_find_out_the_differences = '''
        self.Layout()
        self.Fit()
        self.Refresh()
        self.Update()
        '''


class MyFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, parent=None, title="Add / Remove Buttons")
        self.fSizer = wx.BoxSizer(wx.VERTICAL)
        panel = MyPanel(self)

        self.fSizer.Add(panel, 1, wx.EXPAND)
        self.SetSizer(self.fSizer)
        self.Fit()
        self.Show()

if __name__ == "__main__":
    app = wx.App(False)
    frame = MyFrame()
    app.MainLoop()
