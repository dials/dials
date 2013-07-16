from __future__ import division
from dials.scratch.luiso_s import tst_01, write_2d
from scitbx.array_family import flex
from dials.algorithms.peak_finding import model_2d
import numpy


data2d = numpy.zeros((25, 25), dtype = numpy.float64)
nrow = 6
ncol = 6

ref_ang = 1.0

ref2d = model_2d(nrow, ncol, 1, 2, ref_ang, 55, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[0:6, 0:6] += numpy.float64(data2d_tmp)

ref2d = model_2d(nrow, ncol, 2, 1, ref_ang, 55, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[15:21, 5:11] += numpy.float64(data2d_tmp)

ref2d = model_2d(nrow, ncol, 1, 2, ref_ang, 55, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[18:24, 18:24] += numpy.float64(data2d_tmp)

ref2d = model_2d(nrow, ncol, 2, 1, ref_ang, 55, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[5:11, 17:23] += numpy.float64(data2d_tmp)


ref2d = model_2d(25, 25, 5, 7, ref_ang, 55, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[:, :] += numpy.float64(data2d_tmp)


from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()



'''
data2d = ref2d.as_numpy_array()
'''


print "data2d ="
print data2d

print "________________________________________________________________________"

n = write_2d(flex.double(data2d))
'''
print "a ="
print a

b = numpy.int64(a)
print "b ="
print b
'''
