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

data2d = numpy.zeros((3, 3), dtype = numpy.float64)
data2d[:, :] = 15
data2d[1:2, 1:2] = 50

for row in range(3):
  for col in range(3):
    data2d[row, col] += row * 2
    data2d[row, col] += col * 2

print data2d

mask2d = numpy.zeros((3, 3), dtype = numpy.int32)
mask2d[:, :] = 3
mask2d[1:2, 1:2] = 5

print mask2d

background2d = numpy.copy(data2d)
background2d[:, :] = 0.0

data3d = data2d
data3d.shape = (1,) + data2d.shape
print data3d.shape

mask3d = mask2d
mask3d.shape = (1,) + mask2d.shape
print mask3d.shape

background3d = background2d
background3d.shape = (1,) + background2d.shape
print background3d.shape


from dials.model.data import Reflection, ReflectionList
from scitbx.array_family import flex
r = Reflection()
r.shoebox = flex.double(data2d)
r.shoebox_mask = flex.int(mask2d)
r.shoebox_background = flex.double(background2d)

rlist = ReflectionList()
rlist.append(r)
'''
#from dials.algorithms.background.flat_background_subtractor \
# import tmp_numpy_layering_n_bkgr_avg, layering_and_background_avg
#layering_and_background_avg(rlist)


from dials.algorithms.background.curved_background_subtractor \
layering_and_background_modl
layering_and_background_modl(rlist)
'''

#from dials.algorithms.background.inclined_background_subtractor \
# import layering_and_background_plane
#layering_and_background_plane(rlist)


#tmp_of = '''
from dials.algorithms.integration.summation2d \
 import  flex_2d_layering_n_integrating
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
#'''
