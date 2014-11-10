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

  def __call__(self, np_img_2d, Intst_max, ofst, xyz, title = 'Aaaaaaaaaa',
               mask_colour = False):
    lc_fig = plt.figure()

    #print "imax = ", Intst_max

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
        plt.suptitle(title, fontsize = 18)
        #plt.imshow(np.transpose(np_img_2d), interpolation = "nearest", vmin = 0
        #           , vmax = Intst_max)
        if(mask_colour == False):
          plt.imshow(np.transpose(np_img_2d), interpolation = "nearest",
                     cmap = 'hot', vmin = 0, vmax = Intst_max)
        else:
          plt.imshow(np.transpose(np_img_2d), interpolation = "nearest",
                     cmap = 'gist_rainbow', vmin = 0, vmax = Intst_max)

        xmax = len(np_img_2d[:,1])
        ymax = len(np_img_2d[1,:])
        vl_max = np.amax(np_img_2d)
        vl_min = np.amin(np_img_2d)
        d = vl_max - vl_min
        vl_mid_low = vl_min + d / 3.0
        vl_mid_hig = vl_max - d / 3.0

        if( xmax <= 8 and ymax < 15 ):
          for xpos in range(xmax):
            for ypos in range(ymax):
              f_num = float(np_img_2d[xpos,ypos])
              g = float("{0:.2f}".format(f_num))
              txt_dat = str(g)
              if( g < vl_mid_low ):
                clr_chr = 'yellow'
              elif(g > vl_mid_hig):
                clr_chr = 'black'
              else:
                clr_chr = 'blue'
              plt.annotate(txt_dat, xy = (xpos - 0.3, ypos + 0.3) ,
                           xycoords = 'data', color = clr_chr, size = 22.)


      else:
        plt.suptitle(title, fontsize = 22)
        plt.imshow(np.transpose(np_img_2d), interpolation = "nearest", vmin = 0,
                                cmap = 'gist_rainbow', vmax = 10)

      if(xyz != None):
        arr_w = np.shape(np_img_2d)[0]
        arr_h = np.shape(np_img_2d)[1]

        ini_x = ini_y = - 0.5
        end_y = arr_w - 0.5
        end_x = arr_h - 0.5

        coords = [0,0]
        coords[0] = xyz[0] - 0.5
        coords[1] = xyz[1] - 0.5

        x_lin_from = (ini_x + coords[0]) / 2.0
        y_lin_from = (ini_y + coords[1]) / 2.0

        x_lin_to = (coords[0] + end_x) / 2.0
        y_lin_to = (coords[1] + end_y) / 2.0

        plt.vlines(coords[0], y_lin_from, y_lin_to)
        plt.hlines(coords[1], x_lin_from, x_lin_to)

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
    wx_image = wx_image.Scale(NewW, NewH, wx.IMAGE_QUALITY_NORMAL)
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
