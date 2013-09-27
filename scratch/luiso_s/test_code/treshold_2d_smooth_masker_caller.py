from __future__ import division
import numpy
from iotbx.detectors import ImageFactory


'''
f_names = ["/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0015.cbf", \
           "/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0016.cbf", \
           "/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0017.cbf"]

'''
f_names = ["/dls/i04/data/2013/mx8547-36/Javier/R3_NTD/8_2_collection/8_2_M1S2_2_0015.cbf"]

from dxtbx.sweep import SweepFactory
sweep = SweepFactory.sweep(f_names)
data3d_sw = sweep.to_array()
data3d = data3d_sw.as_numpy_array()
n_frm = numpy.size(data3d[:, 0:1, 0:1])
n_row = numpy.size(data3d[0:1, :, 0:1])
n_col = numpy.size(data3d[0:1, 0:1, :])

data2d = numpy.zeros(n_row * n_col, dtype = numpy.float64).reshape(n_row, n_col)

data2d[:, :] = data3d[0:1, :, :]


from dials.algorithms.peak_finding.spot_finder_lui import do_all_2d


times = 5
shift = 4
n_blocks_x = 12
n_blocks_y = 10
dimensions = "2d"

do_all_2d(sweep, times, shift, n_blocks_x, n_blocks_y, dimensions)


'''

print rlist[0]


from matplotlib import pyplot as plt
print "Plotting data2d"
plt.imshow(data2d, interpolation = "nearest")
plt.show()
'''
