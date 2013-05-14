import numpy
from dials.algorithms.background.flat_subtraction import flat_background_subtraction_2d

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

# def layer_treatment(self, data2d, diffdata2d_ext):
bkgr = flat_background_subtraction_2d(data2d, mask2d)
print 'bkgr =', bkgr
print data2d
