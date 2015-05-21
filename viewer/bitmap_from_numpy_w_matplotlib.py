#!/usr/bin/env python
#
#  bitmap_from_numpy_w_matplotlib.py
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
from dials.algorithms.shoebox import MaskCode

from dials.array_family import flex
from dials_viewer_ext import rgb_img

class wxbmp_from_np_array(object):

  def __init__(self, lst_data_in, show_nums = True, lst_data_mask_in = None):
    if(lst_data_in == [None] and lst_data_mask_in == [None] ):
      self._ini_wx_bmp_lst = None

    else:
      self._ini_wx_bmp_lst = []
      for lst_pos in range(len(lst_data_in)):
        print "lst_pos =", lst_pos
        data_3d_in = lst_data_in[lst_pos]
        xmax = data_3d_in.shape[1]
        ymax = data_3d_in.shape[2]
        # remember to put here some assertion to check that
        # both arrays have the same shape

        if(lst_data_mask_in != None):
          data_3d_in_mask = lst_data_mask_in[lst_pos]

        log_msg = '''
        print "len(data_3d_in) =", len(data_3d_in)
        '''
        self.vl_max = np.amax(data_3d_in)
        self.vl_min = np.amin(data_3d_in)
        tmp_data2d = np.zeros( (xmax, ymax), 'double')
        tmp_data2d_mask = np.zeros( (xmax, ymax), 'double')
        z_dp = data_3d_in.shape[0]
        single_block_lst_01 = []
        for z in range(z_dp):
          print "z =", z
          tmp_data2d[:, :] = data_3d_in[z:z + 1, :, :]
          if(lst_data_mask_in != None):
            tmp_data2d_mask[:, :] = data_3d_in_mask[z:z + 1, :, :]
          else:
            tmp_data2d_mask = None

          #traditional_generator_of_images_with_matplotlib = '''
          data_sigle_img = self._wx_img(tmp_data2d, show_nums, tmp_data2d_mask)
          #'''

          new_generator_of_images_with_cpp = '''
          data_sigle_img = self._wx_img_w_cpp(tmp_data2d, show_nums, tmp_data2d_mask)
          #'''

          single_block_lst_01.append(data_sigle_img)

        self._ini_wx_bmp_lst.append(single_block_lst_01)

  def bmp_lst_scaled(self, scale = 1.0):
    if( self._ini_wx_bmp_lst == None):

      NewW = 350
      NewH = 50

      wx_image = wx.EmptyImage(NewW, NewW)
      wxBitmap = wx_image.ConvertToBitmap()


      dc = wx.MemoryDC(wxBitmap)
      text = 'No Shoebox data'
      w, h = dc.GetSize()
      tw, th = dc.GetTextExtent(text)
      dc.Clear()
      dc.DrawText(text, (w - tw) / 2, (h - th) / 2) #display text in center
      dc.SelectObject(wxBitmap)
      del dc
      wx_bmp_lst = [[wxBitmap]]

    else:
      wx_bmp_lst = []
      for data_3d in self._ini_wx_bmp_lst:
        single_block_lst = []
        for sigle_img_data in data_3d:
          single_block_lst.append(self._wx_bmp_scaled(sigle_img_data, scale))

        wx_bmp_lst.append(single_block_lst)

    return wx_bmp_lst



  def _wx_img_w_cpp(self, np_2d_tmp, show_nums, np_2d_mask = None):
    xmax = np_2d_tmp.shape[0]
    ymax = np_2d_tmp.shape[1]

    if(np_2d_mask == None):
      np_2d_mask = np.zeros( (xmax, ymax), 'double')

    transposed_data = np.zeros( (ymax, xmax), 'double')
    transposed_mask = np.zeros( (ymax, xmax), 'double')

    transposed_data[:,:] = np.transpose(np_2d_tmp)
    transposed_mask[:,:] = np.transpose(np_2d_mask)

    flex_data_in = flex.double(transposed_data)
    flex_mask_in = flex.double(transposed_mask)

    wx_bmp_arr = rgb_img()
    img_array_tmp = wx_bmp_arr.gen_bmp(flex_data_in, flex_mask_in)
    img_array_tmp =    img_array_tmp.as_numpy_array()

    height = np.size( img_array_tmp[:, 0:1, 0:1] )
    width = np.size(  img_array_tmp[0:1, :, 0:1] )
    img_array = np.empty( (height, width, 3),'uint8')
    img_array[:,:,:] = img_array_tmp[:,:,:]

    self._wx_image = wx.EmptyImage(width, height)
    self._wx_image.SetData( img_array.tostring() )

    data_to_become_bmp = (self._wx_image, width, height)

    return data_to_become_bmp



  def _wx_img(self, np_2d_tmp, show_nums, np_2d_mask = None):

    d = self.vl_max - self.vl_min
    vl_mid_low = self.vl_min + d / 3.0
    vl_mid_hig = self.vl_max - d / 3.0

    lc_fig = plt.figure(frameon=False)

    xmax = np_2d_tmp.shape[0]
    ymax = np_2d_tmp.shape[1]

    lc_fig.set_size_inches(xmax * .5, ymax * .5)

    ax = plt.Axes(lc_fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    lc_fig.add_axes(ax)

    a=lc_fig.gca()
    a.set_frame_on(False)
    a.set_xticks([])
    a.set_yticks([])
    plt.axis('off')

    plt.imshow(np.transpose(np_2d_tmp), interpolation = "nearest", cmap = 'hot',
               vmin = self.vl_min, vmax = self.vl_max)

    l_wdt = 1.2

    if( np_2d_mask != None):
      dst_diag = 0.015
      for xpos in range(xmax):
        for ypos in range(ymax):
          loc_mask = int(np_2d_mask[xpos, ypos])

          if( loc_mask != 0 ):

            # drawing box
            plt.vlines(xpos - 0.45, ypos - 0.45, ypos + 0.45,
                       color = 'gray', linewidth=l_wdt)
            plt.vlines(xpos + 0.45, ypos - 0.45, ypos + 0.45,
                       color = 'gray', linewidth=l_wdt)
            plt.hlines(ypos - 0.45, xpos - 0.45, xpos + 0.45,
                       color = 'gray', linewidth=l_wdt)
            plt.hlines(ypos + 0.45, xpos - 0.45, xpos + 0.45,
                       color = 'gray', linewidth=l_wdt)

            if( (loc_mask & MaskCode.Background) == MaskCode.Background ):
              # drawing v_lines
              plt.vlines(xpos, ypos - 0.45, ypos + 0.45,
                         color = 'gray', linewidth=l_wdt)
              plt.vlines(xpos + 0.225, ypos - 0.45, ypos + 0.45,
                         color = 'gray', linewidth=l_wdt)
              plt.vlines(xpos - 0.225, ypos - 0.45, ypos + 0.45,
                         color = 'gray', linewidth=l_wdt)

            if( (loc_mask & MaskCode.BackgroundUsed) == MaskCode.BackgroundUsed ):
              # drawing h_lines
              plt.hlines(ypos , xpos - 0.45, xpos + 0.45,
                         color = 'gray', linewidth=l_wdt)
              plt.hlines(ypos + 0.225, xpos - 0.45, xpos + 0.45,
                         color = 'gray', linewidth=l_wdt)
              plt.hlines(ypos - 0.225, xpos - 0.45, xpos + 0.45,
                         color = 'gray', linewidth=l_wdt)

            if( (loc_mask & MaskCode.Foreground) == MaskCode.Foreground ):
              # drawing lines from 1p5 to 7p5
              plt.plot([xpos - dst_diag, xpos - 0.45], [ypos - 0.45, ypos - dst_diag],
                       color='gray', linewidth=l_wdt)
              plt.plot([xpos + dst_diag, xpos + 0.45], [ypos + 0.45, ypos + dst_diag],
                       color='gray', linewidth=l_wdt)
              plt.plot([xpos - 0.45, xpos + 0.45], [ypos + 0.45, ypos - 0.45],
                       color='gray', linewidth=l_wdt)

            if( (loc_mask &  MaskCode.Valid) == MaskCode.Valid ):
              # drawing lines from 4p5 to 10p5
              plt.plot([xpos - dst_diag, xpos - 0.45], [ypos + 0.45, ypos + dst_diag],
                       color='gray', linewidth=l_wdt)
              plt.plot([xpos + dst_diag, xpos + 0.45], [ypos - 0.45, ypos - dst_diag],
                       color='gray', linewidth=l_wdt)
              plt.plot([xpos + 0.45, xpos - 0.45], [ypos + 0.45, ypos - 0.45],
                       color='gray', linewidth=l_wdt)

    if( show_nums != 0 ):
      for xpos in range(xmax):
        for ypos in range(ymax):
          f_num = np_2d_tmp[xpos,ypos]
          g = float("{0:.3f}".format(float(f_num)))

          txt_dat = str(g)
          if( f_num < vl_mid_low ):
            clr_chr = 'yellow'
          elif(f_num > vl_mid_hig):
            clr_chr = 'black'
          else:
            clr_chr = 'blue'

          plt.annotate(txt_dat, xy = (xpos - 0.5, ypos + 0.1), xycoords = 'data',
                       color = clr_chr, size = 9.)

    lc_fig.canvas.draw()
    width, height = lc_fig.canvas.get_width_height()
    np_buf = np.fromstring (lc_fig.canvas.tostring_rgb(), dtype=np.uint8)
    np_buf.shape = (width, height, 3)
    np_buf = np.roll(np_buf, 3, axis = 2)
    self._wx_image = wx.EmptyImage(width, height)
    self._wx_image.SetData( np_buf )
    data_to_become_bmp = (self._wx_image, width, height)

    plt.close(lc_fig)

    return data_to_become_bmp


  def _wx_bmp_scaled(self, data_to_become_bmp, scale):
    to_become_bmp = data_to_become_bmp[0]
    width = data_to_become_bmp[1]
    height = data_to_become_bmp[2]

    NewW = int(width * scale)
    NewH = int(height * scale)
    to_become_bmp = to_become_bmp.Scale(NewW, NewH, wx.IMAGE_QUALITY_NORMAL)


    wxBitmap = to_become_bmp.ConvertToBitmap()

    imported_from_blog_n_to_consider = '''
    c1 = wx.BitmapButton(self, -1, wx.Bitmap('/home/full/path/to/A.png'))
    '''

    return wxBitmap
