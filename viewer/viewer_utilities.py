#
#  viewer_utilities.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

import numpy as np
import wx
# set backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class np_to_bmp(object):

  def __init__(self):
    print "from init"

  def __call__(self, np_img_2d, Intst_max, ofst, xyz):
    lc_fig = plt.figure()

    if( np_img_2d == None ):

      plt.imshow(np.asarray([[-1]]), interpolation = "nearest")

      lc_fig.canvas.draw()
      self.width, self.height = lc_fig.canvas.get_width_height()
      self.np_buf = np.fromstring ( lc_fig.canvas.tostring_rgb()
                                      , dtype=np.uint8 )
      self.np_buf.shape = (self.width, self.height, 3)
      self.np_buf = np.roll(self.np_buf, 3, axis = 2)
      self.image = wx.EmptyImage(self.width, self.height)
      self.image.SetData( self.np_buf )

      plt.close(lc_fig)
    else:
      if Intst_max > 0:
        plt.imshow(np.transpose(np_img_2d), interpolation = "nearest", vmin = 0
                   , vmax = Intst_max)
        #plt.imshow(np_img_2d, interpolation = "nearest", vmin = 0
        #           , vmax = Intst_max)
      else:
        plt.imshow(np.transpose(np_img_2d), interpolation = "nearest", vmin = 0, vmax = 10)

      if(xyz != None):
        arr_w = np.shape(np_img_2d)[0]
        arr_h = np.shape(np_img_2d)[1]

        #TODO check convention of coordinates DIALS vs Matplotlib

        plt.vlines(xyz[0] - 0.5, xyz[1] / 2.0, (arr_h + xyz[1]) / 2.0)
        plt.hlines(xyz[1] - 0.5, xyz[0] / 2.0, (arr_w + xyz[0]) / 2.0)

      calc_ofst = True
      if(calc_ofst == True):
        ax = lc_fig.add_subplot(1,1,1)

        xlabl = ax.xaxis.get_majorticklocs()
        if(len(xlabl) > 5):
          to_many_labels = True
        else:
          to_many_labels = False
        x_new_labl =[]
        #print xlabl
        for pos in range(len(xlabl)):
          if( float(pos) / 2.0 == int(pos / 2) or to_many_labels == False):
            x_new_labl.append(str(xlabl[pos] + ofst[0] + 0.5))
          else:
            x_new_labl.append("")
        ax.xaxis.set_ticklabels(x_new_labl)

        ylabl = ax.yaxis.get_majorticklocs()
        if(len(ylabl) > 4):
          to_many_labels = True
        else:
          to_many_labels = False
        y_new_labl =[]
        #print ylabl
        for pos in range(len(ylabl)):
          if( float(pos) / 2.0 == int(pos / 2) or to_many_labels == False):
            y_new_labl.append(str(ylabl[pos] + ofst[2] + 0.5))
          else:
            y_new_labl.append("")
        ax.yaxis.set_ticklabels(y_new_labl)

      lc_fig.canvas.draw()
      self.width, self.height = lc_fig.canvas.get_width_height()
      self.np_buf = np.fromstring ( lc_fig.canvas.tostring_rgb()
                                      , dtype=np.uint8 )
      self.np_buf.shape = (self.width, self.height, 3)
      self.np_buf = np.roll(self.np_buf, 3, axis = 2)
      self.image = wx.EmptyImage(self.width, self.height)
      self.image.SetData( self.np_buf )

      plt.close(lc_fig)

    return self.image

  def from_wx_image_to_wx_bitmap(self, wx_image, scale, empty = False):

    NewW = int(self.width * scale)
    NewH = int(self.height * scale)
    wx_image = wx_image.Scale(NewW, NewH, wx.IMAGE_QUALITY_HIGH)
    wxBitmap = wx_image.ConvertToBitmap()

    if(empty == True):
      dc = wx.MemoryDC(wxBitmap)
      text = 'No Image'
      w, h = dc.GetSize()
      tw, th = dc.GetTextExtent(text)
      dc.Clear()
      dc.DrawText(text, (w - tw) / 2, (h - th) / 2) #display text in center
      dc.SelectObject(wxBitmap)
      del dc
    return wxBitmap
