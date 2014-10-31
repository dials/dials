#
#  DIALS viewer
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
from dials.viewer.viewer_frame import ReflectionFrame

class viewer_App(wx.App):
  def OnInit(self):
    self.frame = ReflectionFrame(None, -1, "DIALS Reflections Viewer"
    , wx.DefaultPosition,(550,200))

    self.SetTopWindow(self.frame)
    self.frame.Show(True)
    return True

  def table_in(self, loc_tabl):
    self.frame.tabl_to_frame(loc_tabl)


def next_gen_viewer_test(table):
  print "table(ptr) =", table
  from dials.viewer.next_gen_viewer.multi_3D_slice_viewer import show_3d
  print "Here Yes"


  table_row = table[8037]
  data_flex = table_row['shoebox'].data
  show_3d(data_flex)


  lst_flex_dat = []
  for nm in range(9):
    lst_flex_dat.append(table[nm]['shoebox'].data)
    lst_flex_dat.append(table[nm]['shoebox'].mask)

  show_3d(lst_flex_dat)


if __name__ == "__main__":

  import cPickle as pickle
  import sys
  from dials.array_family import flex
  table = flex.reflection_table.from_pickle(sys.argv[1])

  app = viewer_App(redirect=False)
  app.table_in(table)
  app.MainLoop()

  next_gen_viewer_test(table)
