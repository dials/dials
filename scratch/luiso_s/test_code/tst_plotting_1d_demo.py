from __future__ import division
import numpy
from iotbx.detectors import ImageFactory
from matplotlib import pyplot as plt


f_names = ["/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0015.cbf", \
           "/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0016.cbf", \
           "/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0017.cbf"]


from dxtbx.sweep import SweepFactory
sweep = SweepFactory.sweep(f_names)
data3d_sw = sweep.to_array()
data3d = data3d_sw.as_numpy_array()
n_frm = numpy.size(data3d[:, 0:1, 0:1])
n_row = numpy.size(data3d[0:1, :, 0:1])
n_col = numpy.size(data3d[0:1, 0:1, :])
data2d = numpy.zeros(n_row * n_col, dtype = numpy.float64).reshape(n_row, n_col)

data2d[:, :] = data3d[0:1, :, :]


from dials.algorithms.peak_finding import smooth_2d, smooth_3d
from scitbx.array_family import flex


n_times = 1
data2dsmoth_tmp = smooth_2d(flex.double(data2d), n_times).as_numpy_array()
data2dsmoth = numpy.float64(data2dsmoth_tmp)

n_times = 5
data2dsmoth1_tmp = smooth_2d(flex.double(data2d), n_times).as_numpy_array()
data2dsmoth1 = numpy.float64(data2dsmoth1_tmp)

from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()

plt.imshow(data2dsmoth, interpolation = "nearest")
plt.show()

from matplotlib import pylab

#for pos in range(1585, 1590):
pos = 1585

line = numpy.copy(data2d[:, pos])
line1 = numpy.copy(data2dsmoth[:, pos ])
line2 = numpy.copy(data2dsmoth1[:, pos ])

pylab.subplot2grid((2, 3), (0, 0), colspan = 3)
pylab.plot(line)
pylab.subplot2grid((2, 3), (1, 0), colspan = 3)
pylab.plot(line1)
pylab.show()

pylab.subplot2grid((2, 3), (0, 0), colspan = 3)
pylab.plot(line1)
pylab.subplot2grid((2, 3), (1, 0), colspan = 3)
pylab.plot(line2)
pylab.show()
