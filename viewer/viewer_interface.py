#
#  DIALS viewer_interface
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

def extract_n_show(table):

  from dials.viewer.slice_viewer import show_3d

  lst_nm = range(1, 20)
  flex_dat_frst_lst = []
  flex_dat_seg_lst = []

  for nm in lst_nm:
    # next line might be used later to a test reflection as input data
    # table_row = table[nm]

    flex_dat_frst_lst.append(table[nm]['shoebox'].data)
    flex_dat_seg_lst.append(table[nm]['shoebox'].mask)


  not_needed_for_now = '''
  show_3d(flex_dat_frst_lst[0])
  show_3d(flex_dat_frst_lst[0], flex_dat_seg_lst[0])
  show_3d(flex_dat_frst_lst)
  '''
  show_3d(flex_dat_frst_lst, flex_dat_seg_lst)



if __name__ == "__main__":

  import cPickle as pickle
  import sys
  from dials.array_family import flex
  table = flex.reflection_table.from_pickle(sys.argv[1])

  extract_n_show(table)
