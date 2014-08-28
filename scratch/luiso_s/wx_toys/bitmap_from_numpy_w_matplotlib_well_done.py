#import os, io, Image, PIL
import wx
import numpy as np
# set backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def GetBitmap_from_np_array(data2d):
  lc_fig = plt.figure()


  #lc_fig.set_size_inches(1, 1)
  ax = plt.Axes(lc_fig, [0., 0., 1., 1.])
  ax.set_axis_off()
  lc_fig.add_axes(ax)


  #data2d = np.transpose(data2d_in)
  plt.imshow(np.transpose(data2d), interpolation = "nearest", cmap = 'hot')
  '''
  plt.axis('off')
  #plt.axis('tight')
  #plt.tight_layout(pad=0.0, h_pad=None, w_pad=None)
  plt.tight_layout(pad=0.0)





  '''



  print "len(data2d[:,1]) =", len(data2d[:,1])
  print "len(data2d[1,:]) =", len(data2d[1,:])
  for xpos in range(len(data2d[:,1])):
    for ypos in range(len(data2d[1,:])):
      print "[xpos,ypos] =", [xpos,ypos]
      #txt_dat = str(data2d[xpos,ypos])

      f_num = data2d[xpos,ypos]
      g = float("{0:.2f}".format(f_num))
      txt_dat = str(g)
      plt.annotate(txt_dat, xy = (xpos - 0.3, ypos + 0.3), xycoords = 'data'
                   , color = 'green', size = 12.)

  lc_fig.canvas.draw()
  width, height = lc_fig.canvas.get_width_height()
  np_buf = np.fromstring (lc_fig.canvas.tostring_rgb(), dtype=np.uint8)
  np_buf.shape = (width, height, 3)
  np_buf = np.roll(np_buf, 3, axis = 2)
  wx_image = wx.EmptyImage(width, height)
  wx_image.SetData(np_buf )
  wxBitmap = wx_image.ConvertToBitmap()

  plt.close(lc_fig)

  return wxBitmap

def build_np_img(width=64, height=64):
  data2d = np.zeros( (width, height),'float')
  print "width, height =", width, height
  for x in range(0, width):
    for y in range(0, height):
      data2d[x,y] = np.sqrt(x*x + y*y)
  #data2d[width/4:width*3/4,height/4:height*3/4] = 0
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

    data2d = build_np_img(width=5, height=4)
    print data2d
    bitmap = GetBitmap_from_np_array(data2d)

    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  app = MyApp(redirect=False)
  app.MainLoop()

