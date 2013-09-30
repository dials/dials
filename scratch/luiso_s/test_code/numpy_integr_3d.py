from __future__ import division
#
#  numpy_integr_2d testing code
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
import numpy

data = numpy.zeros((3, 5, 5), dtype = numpy.float64)

data[0, 0, 1] = data[0, 1, 0] = 1.5
data[0, 0, 2] = data[0, 2, 0] = 3.5
data[0, 0, 3] = data[0, 3, 0] = 2.5
data[0, 0, 4] = data[0, 4, 0] = 4.5
data[0, 4, 1] = data[0, 1, 4] = 1.5
data[0, 4, 2] = data[0, 2, 4] = 3.5
data[0, 4, 3] = data[0, 3, 4] = 2.5
data[0, 0, 0] = data[0, 4, 4] = 5.5
data[0, 1:4, 1:4] = 15
data[0, 2:3, 2:3] = 40

data[1, 0, 1] = data[1, 1, 0] = 1.7
data[1, 0, 2] = data[1, 2, 0] = 3.7
data[1, 0, 3] = data[1, 3, 0] = 2.7
data[1, 0, 4] = data[1, 4, 0] = 4.7
data[1, 4, 1] = data[1, 1, 4] = 1.7
data[1, 4, 2] = data[1, 2, 4] = 3.7
data[1, 4, 3] = data[1, 3, 4] = 2.7
data[1, 0, 0] = data[1, 4, 4] = 5.7
data[1, 1:4, 1:4] = 20
data[1, 2:3, 2:3] = 70

data[2, 0, 1] = data[2, 1, 0] = 1.4
data[2, 0, 2] = data[2, 2, 0] = 3.4
data[2, 0, 3] = data[2, 3, 0] = 4.4
data[2, 0, 4] = data[2, 4, 0] = 2.4
data[2, 4, 1] = data[2, 1, 4] = 1.4
data[2, 4, 2] = data[2, 2, 4] = 3.4
data[2, 4, 3] = data[2, 3, 4] = 4.4
data[2, 0, 0] = data[2, 4, 4] = 2.4
data[2, 1:4, 1:4] = 12
data[2, 2:3, 2:3] = 35

print data

mask = numpy.zeros((3, 5, 5), dtype = numpy.int32)
mask[:, :, :] = 3
mask[:, 1:4, 1:4] = 5
print mask

background = numpy.copy(data)
background[:, :, :] = 0.0

from dials.model.data import Reflection, ReflectionList
from scitbx.array_family import flex
r = Reflection()
r.shoebox = flex.double(data)
r.shoebox_mask = flex.int(mask)
r.shoebox_background = flex.double(background)

rlist = ReflectionList()
rlist.append(r)
from dials.algorithms.background.flat_background_subtractor \
 import tmp_numpy_layering_n_bkgr_avg, layering_and_background_avg

from dials.algorithms.background.curved_background_subtractor \
 import tmp_numpy_layering_n_bkgr_modl, layering_and_background_modl

#tmp_numpy_layering_n_bkgr_avg(rlist)
#layering_and_background_avg(rlist)
layering_and_background_modl(rlist)
#tmp_numpy_layering_n_bkgr_modl(rlist)

from dials.algorithms.integration.summation2d \
 import flex_2d_layering_n_integrating

flex_2d_layering_n_integrating(rlist)

for r in rlist:
    print r
    matrix_img = r.shoebox.as_numpy_array()
    print
    print "shoebox"
    print matrix_img

    matrix_bkg = r.shoebox_background.as_numpy_array()
    print
    print "background"
    print matrix_bkg

    matrix_mask = r.shoebox_mask.as_numpy_array()
    print
    print "mask"
    print matrix_mask

