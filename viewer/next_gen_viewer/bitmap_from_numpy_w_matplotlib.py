#!/usr/bin/env python
#
# dials.reflection_viewer.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import wx
import numpy as np
# set backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class wxbmp_from_np_array(object):

  def get_bmp(self, data2d_in):
    print "data2d_in =", data2d_in

    print "data2d_in[0:1, :, 0:1] =", data2d_in[0:1, :, 0:1]
    print "data2d_in[0:1, 0:1, :] =", data2d_in[0:1, 0:1, :]

    xmax = len(data2d_in[0:1, :, 0:1])
    ymax = len(data2d_in[0:1, 0:1, :])

    print data2d.shape()
    print "xmax, ymax =", xmax, ymax

    data2d = np.zeros( (xmax, ymax),'double')
    print "data2d =", data2d
    data2d[:, :] = data2d_in[0:1, :, :]
    print "data2d =", data2d

    vl_max = np.amax(data2d)
    vl_min = np.amin(data2d)
    d = vl_max - vl_min
    vl_mid_low = vl_min + d / 3.0
    vl_mid_hig = vl_max - d / 3.0

    lc_fig = plt.figure(frameon=False)
    lc_fig.set_size_inches(xmax * .6, ymax * .6)
    ax = plt.Axes(lc_fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    lc_fig.add_axes(ax)
    plt.imshow(np.transpose(data2d), interpolation = "nearest", cmap = 'hot')

    print "xmax =", xmax
    print "ymax =", ymax
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
        plt.annotate(txt_dat, xy = (xpos - 0.3, ypos + 0.3), xycoords = 'data'
                     , color = clr_chr, size = 12.)

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
