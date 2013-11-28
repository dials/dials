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

size_box=25
data2d = numpy.zeros((size_box, size_box), dtype = numpy.float64)

data2d[0, 1:5] = data2d[1:5, 0] = 0
data2d[4, 1:4] = data2d[1:4, 4] = 0
data2d[0, 0] = data2d[4, 4] = 0

data2d[6:8, 6:8] = 20
data2d[7, 7] = 60

mod_ncol = 25
mod_nrow = 25
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d
ref2d = model_2d(mod_nrow, mod_ncol, 3, 6, 0.2, 85, 0.5)
data2d_tmp = ref2d.as_numpy_array()

import random


mask2d = numpy.zeros((size_box, size_box), dtype = numpy.int32)
for row in range(size_box):
  for col in range(size_box):
    planepos = row * 0.3 + col * 0.3
    if data2d_tmp[row, col] > planepos:
      data2d[row, col] = data2d_tmp[row, col]
      mask2d[row, col] = 5
    else:
      data2d[row, col] = planepos
      mask2d[row, col] = 3
    data2d[row, col] += random.random() * 6.0
#print data2d
#print mask2d

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
r.shoebox = flex.double(data3d)
r.shoebox_mask = flex.int(mask3d)
r.shoebox_background = flex.double(background3d)

rlist = ReflectionList()
rlist.append(r)
from dials.algorithms.background.flat_background_subtractor \
 import tmp_numpy_layering_n_bkgr_avg, layering_and_background_avg
#layering_and_background_avg(rlist)
#tmp_numpy_layering_n_bkgr_avg(rlist)

from dials.algorithms.background.curved_background_subtractor \
 import tmp_numpy_layering_n_bkgr_modl, layering_and_background_modl \
#tmp_numpy_layering_n_bkgr_modl(rlist)
#layering_and_background_modl(rlist)


from dials.algorithms.background.inclined_background_subtractor \
 import layering_and_background_plane

layering_and_background_plane(rlist)


from dials.algorithms.integration.summation2d \
 import  flex_2d_layering_n_integrating

flex_2d_layering_n_integrating(rlist)




for r in rlist:
  print r
  matrix_img = r.shoebox.as_numpy_array()
  matrix_bkg = r.shoebox_background.as_numpy_array()
  matrix_mask = r.shoebox_mask.as_numpy_array()



from matplotlib import pylab

data2d_tmp =  numpy.copy(matrix_img[0 , :, :])
#print "data2d_tmp ", data2d_tmp
data1d = numpy.copy(data2d_tmp[12, :])
pylab.plot(data1d)

bkg_data2d_tmp =  numpy.copy(matrix_bkg[0 , :, :])
#print "data2d_tmp ", bkg_data2d_tmp
bkg_data1d = numpy.copy(bkg_data2d_tmp[12, :])
pylab.plot(bkg_data1d)
pylab.show()

from matplotlib import pyplot as plt
ind = numpy.arange(25)
p1 = plt.bar(ind, data1d,   0.85, color='r')
plt.show()


print "Plotting shoebox"
plt.imshow(data2d_tmp, interpolation = "nearest")
plt.show()

plt.imshow(bkg_data2d_tmp, interpolation = "nearest")
plt.show()