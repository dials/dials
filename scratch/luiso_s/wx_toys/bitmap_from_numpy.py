import os
import wx
import numpy


def GetBitmap(width=32, height=32, colour = (0,0,0) ):
  array = numpy.zeros( (height, width, 3),'uint8')
  array[:,:,] = colour
  image = wx.EmptyImage(width,height)
  image.SetData( array.tostring())
  wxBitmap = image.ConvertToBitmap()       # OR:  wx.BitmapFromImage(image)
  return wxBitmap

def GetBitmap_from_np_array(data2d):
  height = numpy.size( data2d[0:1, :] )
  width = numpy.size( data2d[:, 0:1] )

  img_array = numpy.zeros( (width, height, 3),'uint8')
  print data2d.max()

  div_scale = 1.0 / data2d.max()
  data2d_scale = numpy.multiply(data2d, div_scale * 255)
  print "div_scale =", div_scale
  print "data2d_scale.max = ", data2d_scale.max()
  data2d_scale[:,:] = data2d_scale.astype(numpy.uint8)
  #a.astype(numpy.int64)
  img_array[:,:,0] = data2d_scale[:,:]
  img_array[:,:,1] = data2d_scale[:,:]
  img_array[:,:,2] = data2d_scale[:,:]

  print img_array.max()
  image = wx.EmptyImage(width,height)
  image.SetData( img_array.tostring())
  wxBitmap = image.ConvertToBitmap()       # OR:  wx.BitmapFromImage(image)
  return wxBitmap
def build_np_img(width=64, height=64):
  data2d = numpy.zeros( (width, height),'double')
  for col in range(0, width):
    for row in range(0, height):
      data2d[col,row] = col * 2 + row * 2
  print data2d.max()
  return data2d

class MyApp(wx.App):
  def OnInit(self):
    self.frame = MyFrame(None, title="Bitmaps")
    self.SetTopWindow(self.frame)
    self.frame.Show()
    return True

class MyFrame(wx.Frame):
  def __init__(self, parent, id=wx.ID_ANY, title="",
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style=wx.DEFAULT_FRAME_STYLE,
               name="MyFrame"):
    super(MyFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)
    # Attributes
    self.panel = wx.Panel(self)

    #bitmap = GetBitmap(width=320, height=50, colour = (100,0,0) )

    data2d = build_np_img(width=300, height=100)
    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(False)
  app.MainLoop()