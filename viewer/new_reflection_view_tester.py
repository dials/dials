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
from dials.array_family import flex

table = flex.reflection_table.from_pickle("all_refl.pickle")

from dials.viewer.tools import show_reflection
show_reflection(table[len(table)//2])
