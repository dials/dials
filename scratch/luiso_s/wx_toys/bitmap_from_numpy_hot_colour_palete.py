#import os
import wx
import numpy


def GetBitmap_from_np_array(data2d):
  width = numpy.size( data2d[0:1, :] )
  height = numpy.size( data2d[:, 0:1] )

  img_array = numpy.zeros( (height ,width, 3),'uint8')
  print "data2d.max =", data2d.max()

  div_scale = 256.0 * 3.0 / data2d.max()
  data2d_scale = numpy.multiply(data2d, div_scale)
  print "div_scale =", div_scale
  print "data2d_scale.max = ", data2d_scale.max()
  #int_data2d_scale = data2d_scale.astype(numpy.uint8)

  for col in range(0, width):
    for row in range(0, height):
      if(data2d_scale[row,col] <= 255):
        img_array[row,col,0] = data2d_scale[row,col]
        img_array[row,col,1] = 0
        img_array[row,col,2] = 0

      elif( data2d_scale[row,col] > 255 and
            data2d_scale[row,col] < 256 * 2):
        img_array[row,col,0] = 255
        img_array[row,col,1] = data2d_scale[row,col] - 256
        img_array[row,col,2] = 0

      else:
        img_array[row,col,0] = 255
        img_array[row,col,1] = 255
        img_array[row,col,2] = data2d_scale[row,col] - 256 * 2 - 1


  print "img_array.max =", img_array.max()
  image = wx.EmptyImage(width,height)
  image.SetData( img_array.tostring())
  wxBitmap = image.ConvertToBitmap()       # OR:  wx.BitmapFromImage(image)
  return wxBitmap
def build_np_img(width=64, height=64):
  data2d = numpy.zeros( (width, height),'float')
  print "width, height =", width, height
  for x in range(0, width):
    for y in range(0, height):
      data2d[x,y] = numpy.sqrt(x*x + y*y)

  data2d[width*1/2:width*3/4,height/4:height*3/4] = data2d.max()
  data2d[width/4:width*1/2,height/4:height*3/4] = 0
  return data2d

class MyApp(wx.App):
  def OnInit(self):
    self.frame = MyFrame(None, title="Bitmaps")
    self.SetTopWindow(self.frame)
    self.frame.Show()
    return True

class MyFrame(wx.Frame):
  def __init__(self, parent, id = wx.ID_ANY, title = "",
               pos = wx.DefaultPosition, size=wx.DefaultSize,
               style = wx.DEFAULT_FRAME_STYLE,
               name="MyFrame"):
    super(MyFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)
    # Attributes
    self.panel = wx.Panel(self)

    data2d = build_np_img(width = 500, height = 1200)
    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(redirect=False)
  app.MainLoop()
