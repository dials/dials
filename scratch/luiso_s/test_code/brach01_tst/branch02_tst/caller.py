from __future__ import division

from dials.scratch.luiso_s.tree_folder_call_test.function01 import *
#from dials.scratch.luiso_s.tree_folder_call_test.function01 import find_mask_2d
import numpy
import time
#import os
from dxtbx.sweep import SweepFactory
import libtbx.load_env

# Try to find the dials regression directory
try:
    dials_regression = libtbx.env.dist_path('dials_regression')
except KeyError, e:
    print 'FAIL: dials_regression not configured'
    quit()



filenames = ["/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0001.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0002.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0003.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0004.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0005.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0006.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0007.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0008.cbf", \
             "/scratch/my_code/dials/dials_data/data_from__dls__i02__data__2013__nt5964-1__2013_02_08_GW__DNA__P1__X4/X4_lots_M1S4_1_0009.cbf"]

sweep = SweepFactory.sweep(filenames)
print "OK 01"

array_01 = sweep.to_array()
data3d = array_01.as_numpy_array()

n_frm = numpy.size(data3d[:, 0:1, 0:1])
n_row = numpy.size(data3d[0:1, :, 0:1])
n_col = numpy.size(data3d[0:1, 0:1, :])
#print 'n_frm, n_col,n_row =', n_frm, n_col, n_row
#print data3d[:, 0:10, 0:10]
#print "OK 02"

data2d = data3d[1, 1000:1500, 1000:1500]
#data2d = data3d[1, :, :]

#display_image_with_predicted_spots_n_centoids(image, xcoords, ycoords, xc, yc):
#"""Display the image with coordinates overlayed."""
from matplotlib import pylab, cm
#plt = pylab.imshow(data2d, cmap = cm.Greys_r, interpolation = 'nearest', origin = 'lower')
#pylab.show()

sumdat = numpy.sum(data2d)
print 'sum =', sumdat

print data2d[206:213, 335:343]



time1 = time.time()
print "time1 =", time1

dif = find_mask_2d(data2d)
lst_box_pos = find_bound_2d(dif)

time2 = time.time()
print "time2 =", time2
timedif = time2 - time1
print "timedif =", timedif

print lst_box_pos

print data2d[206:213, 335:343]
print dif[206:213, 335:343]
sumdat = numpy.sum(data2d)
print 'sum =', sumdat

plt = pylab.imshow(dif, cmap = cm.Greys_r, interpolation = 'nearest', origin = 'lower')
pylab.show()


#pylab.scatter(xcoords, ycoords, marker = 'x')
#pylab.scatter(xc, yc, marker = 'x')
#plt.axes.get_xaxis().set_ticks([])
#plt.axes.get_yaxis().set_ticks([])

