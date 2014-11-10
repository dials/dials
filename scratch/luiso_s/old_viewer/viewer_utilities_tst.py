#
#  viewer_utilities.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

import numpy, wx
# set backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def from_wx_image_to_wx_bitmap(wx_image, width, height, scale):

  NewW = int(width * scale)
  NewH = int(height * scale)
  wx_image = wx_image.Scale(NewW, NewH, wx.IMAGE_QUALITY_HIGH)
  #image.SetData( np_buf.tostring()) # looks like there is no need to convert
  wxBitmap = wx_image.ConvertToBitmap()
  return wxBitmap

  '''
  NewW = int(width * scale)
  NewH = int(height * scale)
  wx_image = wx_image.Scale(NewW, NewH, wx.IMAGE_QUALITY_HIGH)

  wxBitmap = wx_image.ConvertToBitmap()

  dc = wx.MemoryDC(wxBitmap)
  text = 'whatever'
  w, h = dc.GetSize()
  tw, th = dc.GetTextExtent(text)
  dc.DrawText(text, (w - tw) / 2, (h - th) / 2) #display text in center
  dc.SelectObject(wxBitmap)
  del dc
  '''
  return wxBitmap



def build_np_img(width = 64, height = 64):
  '''
  a = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,7,0,7,7,0,0,0,0,7,7,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,9,4,0,0,9,0,0,9,4,0,4,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,9,0,0,0,9,0,0,9,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,9,0,0,0,9,0,0,4,9,7,9,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
       [7,7,7,0,0,0,0,7,7,7,7,7,0,0,0,7,7,7,0,0,0,0,7,7,7,7,0,0,0,7,7,7,0],\
       [0,0,9,0,0,0,0,9,0,9,0,9,0,0,4,0,0,0,9,0,0,9,4,0,4,9,0,0,9,4,0,0,9],\
       [0,0,9,0,0,0,0,9,0,9,0,9,0,0,7,4,4,4,9,0,0,9,0,0,0,9,0,0,9,4,4,4,4],\
       [7,7,9,7,7,0,0,9,0,9,0,9,0,0,4,7,7,4,9,0,0,4,9,7,4,9,0,0,4,9,7,7,4],\
       [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,9,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,0,0,0,0,0,0,0,0],\
       [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
  data2d = numpy.transpose( numpy.asarray(a))
  '''
  data2d = numpy.asarray([[-1]])
  return data2d

class np_to_bmp(object):

  def __init__(self):
    print "from init"

  def __call__(self, np_img_2d, Intst_max, ofst):
    self.fig = plt.figure()
    # remember to make sure this is our convention in (x, y) vs (row, col)


    if Intst_max > 0:
      plt.imshow(numpy.transpose(np_img_2d), interpolation = "nearest", vmin = 0
                 , vmax = Intst_max)
    else:
      plt.imshow(numpy.transpose(np_img_2d), interpolation = "nearest", vmin = 0
                 , vmax = 10)
    calc_ofst = True
    if(calc_ofst == True):
      ax = self.fig.add_subplot(1,1,1)

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

    self.fig.canvas.draw()
    width, height = self.fig.canvas.get_width_height()
    self.np_buf = numpy.fromstring ( self.fig.canvas.tostring_rgb()
                                    , dtype=numpy.uint8 )
    self.np_buf.shape = (width, height, 3)
    self.np_buf = numpy.roll(self.np_buf, 3, axis = 2)
    self.image = wx.EmptyImage(width, height)
    self.image.SetData( self.np_buf )


    plt.close(self.fig)

    return self.image, width, height
