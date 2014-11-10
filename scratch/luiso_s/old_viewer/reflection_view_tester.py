#
#  DIALS viewer test
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

import cPickle as pickle
from dials.viewer.reflection_view import viewer_App
from dials.array_family import flex

table = flex.reflection_table.from_pickle("all_refl.pickle")
#table = flex.reflection_table.from_pickle("/home/lui/cctbx/sources/dials_regression/integration_test_data/i04-weak-data/integrated.pickle")
My_app = viewer_App(redirect=False)
My_app.table_in(table)
My_app.MainLoop()
