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
from __future__ import absolute_import, division, print_function
import wx
from dials.viewer.viewer_low_level_util import (
    flex_arr_img_panel,
    MyGrid,
    flex_3d_frame,
    grid_frame,
)


class show_3d(object):

    """

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

    import six.moves.cPickle as pickle
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

    """

    def __init__(self, flex_arr_one, flex_arr_two=None):
        app = show_3d_wx_app(redirect=False)
        app.in_lst(flex_arr_one, flex_arr_two)
        app.MainLoop()


class show_3d_wx_app(wx.App):
    def OnInit(self):
        self.ImgFrame = flex_3d_frame(None, "3D flex array viewer")
        self.info_panel = flex_arr_img_panel(self.ImgFrame)
        self.ImgFrame.frame_ini_img(self.info_panel)
        return True

    def in_lst(self, flex_lst_one, flex_lst_two=None):
        self.info_panel.ini_n_intro(flex_lst_one, flex_lst_two)
        self.SetTopWindow(self.ImgFrame)
        self.ImgFrame.Show()


class show_reflections(object):
    def __init__(self, table, two_windows=False):

        # two_windows = True
        print("two_windows =", two_windows)

        if two_windows:
            app = show_tabl_2fr_wx_app(redirect=False)
        else:
            app = show_tabl_1fr_wx_app(redirect=False)
        app.in_tabl(table, two_windows)
        app.MainLoop()


class show_tabl_2fr_wx_app(wx.App):
    def OnInit(self):

        self.ImgFrame = flex_3d_frame(None, "DIALS reflections viewer IMG")
        self.flex_panel = flex_arr_img_panel(self.ImgFrame)
        self.ImgFrame.frame_ini_img(self.flex_panel)
        self.ImgFrame.table_exist = True

        self.GridFrame = grid_frame(None, "DIALS reflections viewer Grd")
        self.info_panel = flex_arr_img_panel(self.GridFrame)
        self.data_grid = MyGrid(self.GridFrame)
        self.GridFrame.frame_ini_img(self.flex_panel, self.data_grid)

        return True

    def in_tabl(self, table, two_windows):
        self.data_grid.ini_n_intro(table)
        if two_windows:
            print("two_windows =", two_windows)
            self.flex_panel.ini_n_intro(table)
            self.info_panel.ini_n_intro(table)

            self.GridFrame.Show()
            self.ImgFrame.Show()

        else:
            self.flex_panel.ini_n_intro(table)
            self.ImgFrame.Show()


class show_tabl_1fr_wx_app(wx.App):
    def OnInit(self):
        self.frame = flex_3d_frame(None, "DIALS reflections viewer _")
        self.upper_panel = flex_arr_img_panel(self.frame)
        self.data_grid = MyGrid(self.frame)
        self.frame.frame_ini_img(self.upper_panel, self.data_grid)

        return True

    def in_tabl(self, table, two_windows):

        # if not two_windows:
        self.upper_panel.ini_n_intro(table)
        self.data_grid.ini_n_intro(table)

        self.SetTopWindow(self.frame)
        self.frame.Show()
