import numpy

data2d = numpy.arange(6 * 6 * 6, dtype = numpy.int32).reshape(6, 6, 6)
#mat2d = numpy.arange( ysize * xsize, dtype = 'uintc' ).reshape( ysize, xsize )

data2d[:, 4, 1] += 1
data2d[:, 1, 4] += 1
data2d[:, 4, 2] += 1
data2d[:, 2, 4] += 3
data2d[:, 4, 3] += 1
data2d[:, 3, 4] += 2
data2d[:, 0, 0] += 1
data2d[:, 4, 4] += 5

data2d[:, 1:5, 1:5] += 100
data2d[:, 2:4, 2:4] += 500
print data2d

mask2d = numpy.zeros((6, 6, 6), dtype = numpy.int32)

mask2d[1:5, 1:5, 1:5] = 1

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

from dials.algorithms.integration.integrate2d import Integrate2d

integrate = Integrate2d()

integrate(rlist)
for r in rlist:
    print r
    matrix_img = r.shoebox.as_numpy_array()
    print matrix_img

