from __future__ import division

from dials.scratch.luiso_s.tree_folder_call_test.function01 import *
import numpy
import time
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

array_01 = sweep.to_array()
data3d_ini = array_01.as_numpy_array()
n_frm = numpy.size(data3d_ini[:, 0:1, 0:1])
n_row = numpy.size(data3d_ini[0:1, :, 0:1])
n_col = numpy.size(data3d_ini[0:1, 0:1, :])

print "n_frm,n_row,n_col", n_frm, n_row, n_col


exampl_row_from = 1150
exampl_row_to = 1900
row_dif = exampl_row_to - exampl_row_from


exampl_col_from = 1150
exampl_col_to = 1800
col_dif = exampl_col_to - exampl_col_from

data3d = numpy.zeros((row_dif) * (col_dif) * n_frm , dtype = int).reshape(n_frm, row_dif, col_dif)

data3d = data3d_ini[:, exampl_row_from:exampl_row_to, exampl_col_from:exampl_col_to]

#print 'n_frm, n_col,n_row =', n_frm, n_col, n_row
#print data3d[:, 0:10, 0:10]
#print "OK 02"
#sumdat = numpy.sum(data2d)
#print 'sum =', sumdat
#
#print data2d[206:213, 335:343]
#
#

dif3d = numpy.zeros((row_dif) * (col_dif) * n_frm , dtype = int).reshape(n_frm, row_dif, col_dif)

time1 = time.time()

for frm_tmp in range(n_frm):
    dif3d[frm_tmp, :, :] = find_mask_2d(data3d[frm_tmp, :, :])
time2 = time.time()

dif_3d_ext = find_ext_mask_3d(dif3d)
dif_3d_ext[0:1, :, :] = 0
dif_3d_ext[(n_frm - 1):, :, :] = 0

time3 = time.time()

#for frm_tmp in range(n_frm):
#    data2d = data3d[frm_tmp, :, :
#    x_from_lst, x_to_lst, y_from_lst, y_to_lst = find_bound_2d(dif_3d_ext[frm_tmp, :, :])

x_from_lst, x_to_lst, y_from_lst, y_to_lst, z_from_lst, z_to_lst = find_bound_3d(dif_3d_ext)
time4 = time.time()




for pos in range(len(x_from_lst)):
    print "x_from_lst, x_to_lst, y_from_lst, y_to_lst, z_from_lst, z_to_lst =" \
    , x_from_lst[pos], x_to_lst[pos], y_from_lst[pos], y_to_lst[pos], z_from_lst[pos], z_to_lst[pos]
time5 = time.time()

print "time dif =", time2 - time1
print "time dif =", time3 - time2
print "time dif =", time4 - time3
print "time dif =", time5 - time4

    #print lst_box_pos
    #
    #print data2d[206:213, 335:343]
    #print dif[206:213, 335:343]
    #sumdat = numpy.sum(data2d)
    #print 'sum =', sumdat

from matplotlib import pylab, cm
for frm_tmp in range(n_frm):
    plt = pylab.imshow(dif_3d_ext[frm_tmp, :, :], cmap = cm.Greys_r, interpolation = 'nearest', origin = 'lower')
    pylab.scatter(x_from_lst, y_from_lst, marker = 'x')
    pylab.scatter(x_to_lst, y_to_lst, marker = 'x')
    plt.axes.get_xaxis().set_ticks([])
    plt.axes.get_yaxis().set_ticks([])
    pylab.show()

