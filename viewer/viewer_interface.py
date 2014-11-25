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

  from dials.viewer.slice_viewer import show_reflections

  show_reflections(table)


if __name__ == "__main__":

  import cPickle as pickle
  import sys
  from dials.array_family import flex
  table = flex.reflection_table.from_pickle(sys.argv[1])

  extract_n_show(table)
