from __future__ import division
import boost.python
import numpy as np
import wx
import time

#from looping import loops_2d
from dials.array_family import flex

#from dials_viewer_ext import gen_img, rgb_img
from dials_viewer_ext import rgb_img


def GetBitmap_from_np_array(data2d):


  time1 = time.time()

  #img_array_tmp = gen_img(data2d).as_numpy_array()

  wx_bmp_arr = rgb_img()
  img_array_tmp = wx_bmp_arr.gen_bmp(data2d).as_numpy_array()

  width = np.size(  img_array_tmp[0:1, :, 0:1] )
  height = np.size( img_array_tmp[:, 0:1, 0:1] )
  img_array = np.empty( (height ,width, 3),'uint8')
  img_array[:,:,:] = img_array_tmp[:,:,:]

  time2 = time.time()
  print ("dif(time) =", time2 - time1 )

  print "img_array.max =", img_array.max()
  image = wx.EmptyImage(width,height)
  image.SetData( img_array.tostring())
  wxBitmap = image.ConvertToBitmap()       # OR:  wx.BitmapFromImage(image)
  return wxBitmap


def build_np_img(width=64, height=64):

  data2d = flex.double(flex.grid(width, height), 0)

  for x in range(width):
    for y in range(height):
      data2d[x, y] += (x * 1.5 + y * 1.5)


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

    #data2d = build_np_img(width = 800, height = 1800)
    data2d = build_np_img(width = 142, height = 263)

    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(redirect=False)
  app.MainLoop()
