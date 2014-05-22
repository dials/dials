from __future__ import division
from dials.scratch.luiso_s import hello_tst, tst_01
from scitbx.array_family import flex

x = hello_tst()
print x

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
print "data2d ="
print data2d
a = tst_01(flex.int(data2d)).as_numpy_array()
print "a ="
print a
