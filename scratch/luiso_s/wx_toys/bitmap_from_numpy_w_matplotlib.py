import os, io, Image, PIL
import wx
import numpy as np
from scipy.misc import toimage
import matplotlib.pyplot as plt

#from cStringIO import StringIO

def GetBitmap_from_np_array(data2d):
#def OnExternalUpdate (self, event) :
  """The OnExternalUpdate() function updates the image and the title
  from @p event.
  """

  to_study = '''
  # See self.load_image().
  self._img = event.img
  self.viewer.set_image(self._img)
  if self.settings_frame is not None :
    self.settings_frame.set_image(self._img)
  if self.zoom_frame is not None:
    self.zoom_frame.set_image(self._img)
  self.SetTitle(event.title)
  self.update_statusbar()
  self.Layout()
  if event.spectrum is not None:
    dc=wx.PaintDC(self.panelspectrum)

    img=toimage(data2d,high=256, low=0)

#     im = ImageEnhance.Contrast(img)
#      img =im.enhance(0.3)i
    if self.angle!=0:
      img=self.ScaleRotateTranslate(img, 25, (512,512), None, (-0.6,-0.9),True)
      img=img.rotate(-90)
      img=img.crop((55, 400, 976, 624))
    else: img=img.rotate(90)
    wximg=wx.EmptyImage(img.size[0],img.size[1])
    wximg.SetData(img.convert("RGB").tostring())
    wximg=wximg.Scale(600,300, wx.IMAGE_QUALITY_HIGH)
    bitmp=wx.BitmapFromImage(wximg, depth=-1)
    dc.DrawBitmap(bitmp,50, 0)
  self.draw(event.data,event.number)
  '''



  one_more_way = '''
  from scipy.misc import toimage
  toimage(data2d).show()
  '''
  #Image.fromarray(a)

  #new_way= '''
  np_max = np.amax(data2d)
  print "np_max =", np_max
  img=toimage(data2d,high=256, low=0, pal=3, mode='P')
  #img = Image.fromarray(data2d, 'RGB')
  wximg=wx.EmptyImage(img.size[0],img.size[1])
  wximg.SetData(img.convert("RGB").tostring())
  #wximg.SetData(img.tostring())
  #wximg=wximg.Scale(600,300, wx.IMAGE_QUALITY_HIGH)
  bitmp=wx.BitmapFromImage(wximg, depth=-1)
  wxBitmap = bitmp
  #'''

  old_way = '''
  plt.imshow(data2d, interpolation = "nearest")
  plt.savefig("/dev/shm/img_tmp.png", format = 'png')
  wxBitmap = wx.Bitmap("/dev/shm/img_tmp.png")
  #'''

  return wxBitmap


def build_np_img(width=64, height=64):
  data2d = np.zeros( (width, height),'float')
  print "width, height =", width, height
  for x in range(0, width):
    for y in range(0, height):
      data2d[x,y] = np.sqrt(x*x + y*y)
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

    data2d = build_np_img(width=10, height=14)
    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(redirect=False)
  app.MainLoop()

