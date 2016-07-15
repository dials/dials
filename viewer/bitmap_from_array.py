#!/usr/bin/env python
#
#  bitmap_from_array.py
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

no_longer_needeed_deps = '''
# set backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
'''

from dials.algorithms.shoebox import MaskCode

from dials.array_family import flex
from dials_viewer_ext import rgb_img

class wxbmp_from_np_array(object):

  def __init__(self, lst_data_in, show_nums = True, palette = "black2white", lst_data_mask_in = None):
    self.wx_bmp_arr = rgb_img()
    if lst_data_in == [None] and lst_data_mask_in == [None]:
      self._ini_wx_bmp_lst = None

    else:
      self._ini_wx_bmp_lst = []
      for lst_pos in range(len(lst_data_in)):
        #print "lst_pos =", lst_pos
        data_3d_in = lst_data_in[lst_pos]
        xmax = data_3d_in.shape[1]
        ymax = data_3d_in.shape[2]
        # remember to put here some assertion to check that
        # both arrays have the same shape

        if lst_data_mask_in is not None:
          data_3d_in_mask = lst_data_mask_in[lst_pos]

        log_msg = '''
        print "len(data_3d_in) =", len(data_3d_in)
        '''
        self.vl_max = float(np.amax(data_3d_in))
        self.vl_min = float(np.amin(data_3d_in))
        tmp_data2d = np.zeros((xmax, ymax), 'double')
        tmp_data2d_mask = np.zeros((xmax, ymax), 'double')
        z_dp = data_3d_in.shape[0]
        single_block_lst_01 = []
        for z in range(z_dp):
          #print "z =", z
          tmp_data2d[:, :] = data_3d_in[z:z + 1, :, :]
          if lst_data_mask_in is not None:
            tmp_data2d_mask[:, :] = data_3d_in_mask[z:z + 1, :, :]
          else:
            tmp_data2d_mask = None

          data_sigle_img = self._wx_img_w_cpp(tmp_data2d, show_nums, palette, tmp_data2d_mask)

          single_block_lst_01.append(data_sigle_img)

        self._ini_wx_bmp_lst.append(single_block_lst_01)


  def bmp_lst_scaled(self, scale = 1.0):
    if self._ini_wx_bmp_lst == None:

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


  def _wx_img_w_cpp(self, np_2d_tmp, show_nums, palette, np_2d_mask = None):

    xmax = np_2d_tmp.shape[1]
    ymax = np_2d_tmp.shape[0]

    if np_2d_mask is None:
      np_2d_mask = np.zeros((ymax, xmax), 'double')

    transposed_data = np.zeros((ymax, xmax), 'double')
    transposed_mask = np.zeros((ymax, xmax), 'double')

    transposed_data[:,:] = np_2d_tmp
    transposed_mask[:,:] = np_2d_mask

    flex_data_in = flex.double(transposed_data)
    flex_mask_in = flex.double(transposed_mask)

    err_code = self.wx_bmp_arr.set_min_max(self.vl_min, self.vl_max)

    test_log = '''
    print "self.vl_min, self.vl_max = ", self.vl_min, self.vl_max
    print "err_code =", err_code
    print "before crash"
    print "palette =", palette
    '''

    if palette == "black2white":
      palette_num = 1
    elif palette == "white2black":
      palette_num = 2
    elif palette == "hot ascend":
      palette_num = 3
    else: # assuming "hot descend"
      palette_num = 4

    img_array_tmp = self.wx_bmp_arr.gen_bmp(flex_data_in, flex_mask_in, show_nums, palette_num)

    np_img_array = img_array_tmp.as_numpy_array()

    height = np.size(np_img_array[:, 0:1, 0:1])
    width = np.size( np_img_array[0:1, :, 0:1])
    img_array = np.empty((height, width, 3),'uint8')
    img_array[:,:,:] = np_img_array[:,:,:]

    self._wx_image = wx.EmptyImage(width, height)
    self._wx_image.SetData(img_array.tostring())

    data_to_become_bmp = (self._wx_image, width, height)

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
