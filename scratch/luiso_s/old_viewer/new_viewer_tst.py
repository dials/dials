#
#  DIALS viewer_tester
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
import wx


def next_gen_viewer_test(table):

  from dials.viewer.slice_viewer import show_3d
  table_row = table[8037]
  data_flex = table_row['shoebox'].data
  show_3d(data_flex)


  lst_flex_dat = []
  for nm in range(312000, 312020):
    lst_flex_dat.append(table[nm]['shoebox'].data)
    lst_flex_dat.append(table[nm]['shoebox'].mask)

  show_3d(lst_flex_dat)


if __name__ == "__main__":

  import sys
  from dials.array_family import flex
  from reflection_view import viewer_App

  table = flex.reflection_table.from_pickle(sys.argv[1])

  app = viewer_App()
  app.table_in(table)
  app.MainLoop()

  next_gen_viewer_test(table)
