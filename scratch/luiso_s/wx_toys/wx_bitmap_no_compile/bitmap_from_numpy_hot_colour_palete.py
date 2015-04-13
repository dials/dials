import wx
import numpy as np
import time

from looping import loops_2d

def GetBitmap_from_np_array(data2d):
  width = np.size( data2d[0:1, :] )
  height = np.size( data2d[:, 0:1] )


  print "data2d.max =", data2d.max()

  div_scale = 764.0 / data2d.max()
  data2d_scale = np.multiply(data2d, div_scale)
  print "div_scale =", div_scale
  print "data2d_scale.max = ", data2d_scale.max()

  img_array = np.empty( (height ,width, 3),'uint8')

  time1 = time.time()

  img_array_tmp = loops_2d(data2d_scale)
  img_array[:,:,:] = img_array_tmp[:,:,:]

  time2 = time.time()
  print ("dif(time) =", time2 - time1 )

  print "img_array.max =", img_array.max()
  image = wx.EmptyImage(width,height)
  image.SetData( img_array.tostring())
  wxBitmap = image.ConvertToBitmap()       # OR:  wx.BitmapFromImage(image)
  return wxBitmap


def build_np_img(width=64, height=64):
  data2d = np.zeros( (width, height),'float')
  print "width, height =", width, height
  for x in range(0, width):
    for y in range(0, height):
      data2d[x,y] = np.sqrt(x*x + y*y)

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
               pos = wx.DefaultPosition, size=(1900,900),
               style = wx.DEFAULT_FRAME_STYLE,
               name="MyFrame"):
    super(MyFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)
    # Attributes
    self.panel = wx.Panel(self)

    data2d = build_np_img(width = 800, height = 1800)
    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(redirect=False)
  app.MainLoop()
