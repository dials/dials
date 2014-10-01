#
#  DIALS Image Utilities
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."
#
#
from __future__ import division
import wx
import wx.lib.scrolledpanel as scroll_pan
import numpy as np
# set backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def GetBitmap_from_np_array(data2d):
  lc_fig = plt.figure(frameon=False)

  xmax = len(data2d[:,1])
  ymax = len(data2d[1,:])
  vl_max = np.amax(data2d)
  vl_min = np.amin(data2d)
  d = vl_max - vl_min
  vl_mid_low = vl_min + d / 3.0
  vl_mid_hig = vl_max - d / 3.0

  lc_fig.set_size_inches(xmax * .6, ymax * .6)
  #lc_fig.set_size_inches(xmax * .2, ymax * .2)

  ax = plt.Axes(lc_fig, [0., 0., 1., 1.])
  ax.set_axis_off()
  lc_fig.add_axes(ax)
  plt.imshow(np.transpose(data2d), interpolation = "nearest", cmap = 'hot')

  #print "xmax =", xmax
  #print "ymax =", ymax
  for xpos in range(xmax):
    for ypos in range(ymax):
      f_num = data2d[xpos,ypos]
      g = float("{0:.2f}".format(float(f_num)))

      txt_dat = str(g)
      if( g < vl_mid_low ):
        clr_chr = 'yellow'
      elif(g > vl_mid_hig):
        clr_chr = 'black'
      else:
        clr_chr = 'blue'
      #plt.annotate(txt_dat, xy = (xpos - 0.3, ypos + 0.3), xycoords = 'data'
      #             , color = clr_chr, size = 12.)

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

class ImageListCtrl(scroll_pan.ScrolledPanel):
  """Simple control to display a list of images"""
  def __init__(self, parent, orient, bitmaps=list(),
               style=wx.TAB_TRAVERSAL|wx.BORDER_SUNKEN):
    super(ImageListCtrl, self).__init__(parent,
                                        style = style)

    # Attributes
    self.images = list()
    self.orient = orient.lower()
    if( self.orient  == "portrait" ):
      self.sizer = wx.BoxSizer(wx.VERTICAL)
    else:
      self.sizer = wx.BoxSizer(wx.HORIZONTAL)

    # Setup
    for bmp in bitmaps:
        self.AppendBitmap(bmp)
    self.SetSizer(self.sizer)

  def AppendBitmap(self, bmp):
    """Add another bitmap to the control"""
    self.images.append(bmp)
    sbmp = wx.StaticBitmap(self, bitmap = bmp)
    #self.sizer.Add(sbmp, 0, wx.EXPAND|wx.TOP, 5)
    if( self.orient  == "portrait" ):
      self.sizer.Add(sbmp, 5, wx.EXPAND|wx.TOP, 5)
    else:
      self.sizer.Add(sbmp, 5, wx.EXPAND|wx.LEFT, 5)

    self.SetupScrolling()
