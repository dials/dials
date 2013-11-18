from __future__ import division
import numpy

data2d = numpy.zeros((5, 5), dtype = numpy.int32)
data2d[0, 1] = data2d[1, 0] = 1
data2d[0, 2] = data2d[2, 0] = 3
data2d[0, 3] = data2d[3, 0] = 2
data2d[0, 4] = data2d[4, 0] = 4
data2d[4, 1] = data2d[1, 4] = 1
data2d[4, 2] = data2d[2, 4] = 3
data2d[4, 3] = data2d[3, 4] = 2
data2d[0, 0] = data2d[4, 4] = 5
data2d[1:4, 1:4] = 10
data2d[2:3, 2:3] = 50
print data2d

mask2d = numpy.zeros((5, 5), dtype = numpy.int32)
mask2d[1:4, 1:4] = 1
print mask2d


data3d = data2d
data3d.shape = (1,) + data2d.shape
print data3d.shape

mask3d = mask2d
mask3d.shape = (1,) + mask2d.shape
print mask3d.shape


from dials.model.data import Reflection, ReflectionList
from scitbx.array_family import flex
r = Reflection()
r.shoebox = flex.int(data2d)
r.shoebox_mask = flex.int(mask2d)

rlist = ReflectionList()
rlist.append(r)

from dials.algorithms.centroid.toy_centroid_Lui import toy_centroid_lui

toy_centroid_lui(rlist)

for r in rlist:
  print r
  matrix_img = r.shoebox.as_numpy_array()
  print matrix_img
print "done"
