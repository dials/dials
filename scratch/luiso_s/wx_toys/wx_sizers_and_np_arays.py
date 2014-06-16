#
#  DIALS viewer
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

import wx, os, numpy
import matplotlib.pyplot as plt
class TestFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)

        self.MaxImageSize = 300

        btn_nxt_refl = wx.Button(self, -1, "Next Reflection ")
        b_a = wx.Button(self, -1, "Previous Reflection")
        btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext)
        b_a.Bind(wx.EVT_BUTTON, self.DisplayPrev)

        # starting with an EmptyBitmap
        self.Image = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                     self.MaxImageSize, self.MaxImageSize))

        self.DisplayNext()

        # Using a Sizer to handle the layout: is not recommended to use absolute # positioning

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,10)

        # adding stretchable space before and after centers the image.
        box.Add((1,1),1)
        box.Add(self.Image
                , 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE
                , 10)

        box.Add((1,1),1)
        box.Add(b_a, 0, wx.CENTER | wx.ALL,10)

        self.SetSizerAndFit(box)

        wx.EVT_CLOSE(self, self.OnCloseWindow)

    def DisplayNext(self, event = None):
        print "test 01"
        np_img = build_np_img(width = 64, height = 164)

        Img = GetBitmap_from_np_array(np_img)

        self.Image.SetBitmap(Img)
        self.Fit()
        self.Layout()
        self.Refresh()

    def DisplayPrev(self, event = None):
        print "test 01"
        np_img = build_np_img(width = 164, height = 64)

        Img = GetBitmap_from_np_array(np_img)

        self.Image.SetBitmap(Img)
        self.Fit()
        self.Layout()
        self.Refresh()

    def OnCloseWindow(self, event):
        self.Destroy()

def GetBitmap_from_np_array(np_img_2d):
  fig = plt.figure()
  # remember to make sure this is our convention in (x, y) vs (row, col)
  plt.imshow(numpy.transpose(np_img_2d), interpolation = "nearest")
  fig.canvas.draw()
  width, height = fig.canvas.get_width_height()
  np_buf = numpy.fromstring ( fig.canvas.tostring_rgb(), dtype=numpy.uint8 )
  np_buf.shape = (width, height, 3)
  np_buf = numpy.roll(np_buf, 3, axis = 2)
  image = wx.EmptyImage(width, height)
  image.SetData( np_buf )
  #image.SetData( np_buf.tostring()) # looks like there is no need to convert
  wxBitmap = image.ConvertToBitmap()

  return wxBitmap

def build_np_img(width = 64, height = 64):
  data2d = numpy.zeros( (width, height), 'float')
  print "width, height =", width, height
  for x in range(0, width):
    for y in range(0, height):
      #data2d[x,y] = numpy.sqrt(x*x + y*y)
      data2d[x,y] = x + y
  data2d[width/4:width*3/4,height/4:height*3/4] = 0
  print "data2d.max =", data2d.max()
  return data2d

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