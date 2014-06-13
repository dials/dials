import wx
import numpy
import matplotlib.pyplot as plt

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

    data2d = build_np_img(width = 60, height = 20)
    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(redirect=False)
  app.MainLoop()

