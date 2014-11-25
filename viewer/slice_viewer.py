#!/usr/bin/env python
#
#  slice_viewer.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


from __future__ import division
import wx
from viewer_low_level_util import flex_3d_frame, flex_arr_img_panel


class show_3d(object):


  '''

  This is a useful class for developers to view 3D flex array(s) at low level code

  It opens a WxPython `Frame` with sliced layers of the array(s).

  Depending on the data entered it should show different amount of information

  There are 4 ways to use show_3d():

  1) show_3d(3D_flex_array)
  2) show_3d(3D_flex_array, 3D_flex_array_with_mask)
  3) show_3d(python_list_of_3D_flex_arrays)
  3) show_3d(python_list_of_3D_flex_arrays, python_list_of_3D_flex_arrays_with_mask)

  Basically the idea is:

  Given one single argument it should show shoebox(es)

  or

  Given two arguments it should show shoebox(es) with the mask superimposed, assuming
  the first argument is `shoebox.data` and second argument is shoebox.mask

  Following there is an example of how to see a set of shoebox(es) from a pickle file,
  this is taken from  a test piece of code:

  import cPickle as pickle
  from dials.array_family import flex
  from dials.viewer.slice_viewer import show_3d


  table = flex.reflection_table.from_pickle("PATH/TO/MY/PICKLE/FILE")

  from dials.viewer.slice_viewer import show_3d

  flex_dat_frst_lst = []
  flex_dat_seg_lst = []

  for nm in range(Beginning_from_row , ending_at_row):

    flex_dat_frst_lst.append(table[nm]['shoebox'].data)
    flex_dat_seg_lst.append(table[nm]['shoebox'].mask)

  #to see data from a single shoebox
  show_3d(flex_dat_frst_lst[row])

  #to see data and mask from a single shoebox
  show_3d(flex_dat_frst_lst[row], flex_dat_seg_lst[0])

  #to see data from a set shoeboxes
  show_3d(flex_dat_frst_lst)

  #to see data and mask from a set shoeboxes
  show_3d(flex_dat_frst_lst, flex_dat_seg_lst)

  '''

  def __init__(self, flex_arr_one, flex_arr_two = None):
    app = show_3d_wx_app(redirect=False)
    app.in_lst(flex_arr_one, flex_arr_two)
    app.MainLoop()

class show_reflections(show_3d):
  def __init__(self, table):
    lst_nm = range(1, 20)
    flex_dat_frst_lst = []
    flex_dat_seg_lst = []

    for nm in lst_nm:
      # next line might be used later to a test reflection as input data
      # table_row = table[nm]

      flex_dat_frst_lst.append(table[nm]['shoebox'].data)
      flex_dat_seg_lst.append(table[nm]['shoebox'].mask)


    show_3d(flex_dat_frst_lst, flex_dat_seg_lst)
    # testing log
    #print "table[0] =", table[0]

    output_full_row = '''

    table[0] = {
    'imageset_id': 0,
    'xyzcal.mm': (0.0, 0.0, 0.0),
    'intensity.sum.value': 15.0,
    'xyzobs.px.variance': (0.286144578313253, 0.2676706827309237, 0.25),
    'intensity.sum.variance': 15.0,
    's1': (-0.03497287541155274, 0.604366929910932, -0.826295905178965),
    'shoebox': <dials_model_data_ext.Shoebox object at 0x7b4f820>,
    'iobs': 0,
    'rlp': (-0.035319267434222125, 0.27946863045406983, -0.5712908827798997),
    'xyzobs.mm.variance': (0.00846530120481, 0.007918769477911, 1.7134729863002e-06),
    'miller_index': (36, -4, -17),
    'flags': 0,
    'bbox': (1158, 1161, 68, 70, 0, 1),
    'xyzobs.mm.value': (199.43804512406385, 11.977194484339728, 1.4324789835743459),
    'xyzobs.px.value': (1159.5, 69.23333333333333, 0.5),
    'id': 0,
    'panel': 0
    }

    '''

    not_needed_for_now = '''
    show_3d(flex_dat_frst_lst[0])
    show_3d(flex_dat_frst_lst[0], flex_dat_seg_lst[0])
    show_3d(flex_dat_frst_lst)
    show_3d(flex_dat_frst_lst, flex_dat_seg_lst)
    '''


class show_3d_wx_app(wx.App):
  def OnInit(self):
    self.frame = flex_3d_frame(None, '3D flex array viewer')
    self.upper_panel = flex_arr_img_panel(self.frame)
    self.frame.frame_ini_img(self.upper_panel)
    return True


  def in_lst(self, flex_lst_one, flex_lst_two = None):
    self.upper_panel.ini_n_intro(flex_lst_one, flex_lst_two)
    self.SetTopWindow(self.frame)
    self.frame.Show()

