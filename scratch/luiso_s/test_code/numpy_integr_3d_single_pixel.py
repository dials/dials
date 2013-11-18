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

data = numpy.zeros((4, 4, 4), dtype = numpy.float64)
data[:, :, :] = 5
data[1:3, 1:3, 1:3] = 50
print data
for frame in range(4):
  for row in range(4):
    for col in range(4):
      data[frame, row, col] += frame + col + row
'''
data = numpy.zeros((3, 3, 3), dtype = numpy.float64)
data[:, :, :] = 5
data[1:2, 1:2, 1:2] = 50
print data
for frame in range(3):
    for row in range(3):
        for col in range(3):
            data[frame, row, col] += frame + col + row
'''


mask = numpy.zeros((4, 4, 4), dtype = numpy.int32)
mask[:, :, :] = 3
mask[1:3, 1:3, 1:3] = 5
print mask

background = numpy.copy(data)
background[:, :, :] = 0.0



'''
mask = numpy.zeros((3, 3, 3), dtype = numpy.int32)
mask[:, :, :] = 3
mask[:, 1:2, 1:2] = 5
print mask

background = numpy.copy(data)
background[:, :, :] = 0.0
'''
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

from dials.algorithms.background.inclined_background_subtractor \
 import layering_and_background_plane


#tmp_numpy_layering_n_bkgr_avg(rlist)
#layering_and_background_avg(rlist)
#layering_and_background_modl(rlist)
#tmp_numpy_layering_n_bkgr_modl(rlist)


layering_and_background_plane(rlist)

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
