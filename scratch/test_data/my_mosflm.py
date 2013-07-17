from __future__ import division
from dials.scratch.luiso_s import tst_01, write_2d
from scitbx.array_family import flex
from dials.algorithms.peak_finding import model_2d
from dials.algorithms.peak_finding import smooth_2d, smooth_3d
from matplotlib import pyplot as plt
from matplotlib import pylab

import numpy

data2d = numpy.zeros((50, 80), dtype = numpy.float64)
nrow = 10
ncol = 10



for xpos in range(3):
    for ypos in range(3):
        row_str = ypos * 12
        col_str = xpos * 20
        ref_ang = float(ypos / 10)
        '''flex_int model_2d(int nrow, int ncol, float a, float b,
        float delta_ang, float imax, float asp)'''
        ref2d = model_2d(nrow, ncol, 3, 2, ref_ang, 955, 0.5)
        data2d_tmp = ref2d.as_numpy_array()
        data2d[row_str:row_str + nrow, col_str:col_str + ncol] += numpy.float64(data2d_tmp)

ref2d = model_2d(50, 80, 25, 40, 1.0, 55, 0.5)
data2d_tmp = ref2d.as_numpy_array()
data2d[:, :] += numpy.float64(data2d_tmp)


#n = write_2d(flex.double(data2d))



print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()


'''
n_times = 1
data2dsmoth_tmp = smooth_2d(flex.double(data2d), n_times).as_numpy_array()
data2dsmoth = numpy.float64(data2dsmoth_tmp)


n_times = 5
data2dsmoth_2t_tmp = smooth_2d(flex.double(data2d), n_times).as_numpy_array()
data2dsmoth_2t = numpy.float64(data2dsmoth_2t_tmp)


#n = write_2d(flex.double(data2dsmoth))

print "Plotting data2dsmoth"
plt.imshow(data2dsmoth, interpolation = "nearest")
plt.show()


line = numpy.copy(data2d[230, :])
pylab.subplot2grid((3, 3), (0, 0), colspan = 3)
pylab.plot(line)

line1 = numpy.copy(data2dsmoth[230, :])
pylab.subplot2grid((3, 3), (1, 0), colspan = 3)
pylab.plot(line1)

line2 = numpy.copy(data2dsmoth_2t[230, :])
pylab.subplot2grid((3, 3), (2, 0), colspan = 3)
pylab.plot(line2)

pylab.show()
'''
