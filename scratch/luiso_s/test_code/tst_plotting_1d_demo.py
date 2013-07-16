from __future__ import division
import numpy
from iotbx.detectors import ImageFactory
from matplotlib import pyplot as plt
'''cad = '/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0005.cbf'
/dls/i24/data/2013/mx8423-66/CPV17/./grid1/grid1_1_0001.cbf'''
#cad = '/scratch/my_code/dials/dials_data/xrd_2d_curv_background/grid1_1_0115.cbf'


cad = '/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0015.cbf'
print cad
image = ImageFactory(cad)
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

data = image.get_raw_data()
print 'here 1'

data2d = numpy.zeros(nfast * nslow, dtype = numpy.float64).reshape(nfast, nslow)
for f in range(nfast):
    for s in range(nslow):
        data2d[f, s] = data[s * nfast + f]

print "nslow, nfast =", nslow, nfast



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

from matplotlib import pylab, cm

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
