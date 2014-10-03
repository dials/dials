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

  def get_bmp_lst(self, data_3d_in):
    print "data_3d_in =\n", data_3d_in

    print "data_3d_in[0:1, :, 0:1] =\n", data_3d_in[0:1, :, 0:1]
    print "data_3d_in[0:1, 0:1, :] =\n", data_3d_in[0:1, 0:1, :]

    print "data_3d_in.shape[0] =\n", data_3d_in.shape[0]
    print "data_3d_in.shape[1] =\n", data_3d_in.shape[1]
    print "data_3d_in.shape[2] =\n", data_3d_in.shape[2]

    z_dp = data_3d_in.shape[0]
    self.xmax = data_3d_in.shape[1]
    self.ymax = data_3d_in.shape[2]


    self.vl_max = np.amax(data_3d_in)
    self.vl_min = np.amin(data_3d_in)

    wx_bmp_lst = []
    tmp_data2d = np.zeros( (self.xmax, self.ymax), 'double')
    for z in range(z_dp):
      tmp_data2d[:, :] = data_3d_in[z:z + 1, :, :]
      wx_bmp_lst.append(self.wxbmp(tmp_data2d))

    return wx_bmp_lst

  def wxbmp(self, np_2d_tmp):

    d = self.vl_max - self.vl_min
    vl_mid_low = self.vl_min + d / 3.0
    vl_mid_hig = self.vl_max - d / 3.0
    lc_fig = plt.figure(frameon=False)
    lc_fig.set_size_inches(self.xmax * .6, self.ymax * .6)
    ax = plt.Axes(lc_fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    lc_fig.add_axes(ax)
    plt.imshow(np.transpose(np_2d_tmp), interpolation = "nearest", cmap = 'hot')

    for xpos in range(self.xmax):
      for ypos in range(self.ymax):
        f_num = np_2d_tmp[xpos,ypos]
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
