#
#  DIALS viewer test
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

from dials.viewer.reflection_view import viewer_App
import cPickle as pickle
from dials.array_family import flex

table = flex.reflection_table.from_pickle("all_refl.pickle")
My_app = viewer_App(redirect=False)
My_app.table_in(table)
My_app.MainLoop()
