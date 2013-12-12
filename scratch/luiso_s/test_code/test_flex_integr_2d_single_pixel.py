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
# import numpy
from scitbx.array_family import flex
#data2d = numpy.zeros((3, 3), dtype = numpy.float64)
data2d = flex.double(flex.grid(1, 3, 3),15)
data2d[0, 1, 1] = 50

for row in range(3):
  for col in range(3):
    data2d[0,row, col] += row * 2
    data2d[0,row, col] += col * 2

mask2d = flex.int(flex.grid(1, 3, 3),3)
mask2d[0, 1, 1] = 5

background2d = flex.double(flex.grid(1, 3, 3),0)
from dials.model.data import Reflection, ReflectionList
from scitbx.array_family import flex
r = Reflection()
r.shoebox = data2d
r.shoebox_mask = mask2d
r.shoebox_background = background2d

rlist = ReflectionList()
rlist.append(r)

#from dials.algorithms.background.flat_background_subtractor \
# import layering_and_background_avg
#layering_and_background_avg(rlist)

#from dials.algorithms.background.curved_background_subtractor \
# import layering_and_background_modl
#layering_and_background_modl(rlist)

from dials.algorithms.background.inclined_background_subtractor \
 import layering_and_background_plane
layering_and_background_plane(rlist)


from dials.algorithms.integration.summation2d \
 import  flex_2d_layering_n_integrating
flex_2d_layering_n_integrating(rlist)

for r in rlist:
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

  print r
  if r.intensity == 35:
    print "Summation integration  ...  OK"
