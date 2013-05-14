# toy centroid implementation: this is deliberately bad to illustrate how the
# interface works and how things could be improved
# from __future__ import division
from dials.interfaces.centroid.centroid_interface_prototype import \
    centroid_interface_prototype as centroid_interface

from dials.interfaces.centroid.centroid_interface_prototype import \
    CentroidException

import numpy

class toy_centroid_lui(centroid_interface):
    def __init__(self, reflections):

        centroid_interface.__init__(self, reflections)

        return

    def compute_shoebox_centroid(self, shoebox):
        #import time

        f_size, r_size, c_size = shoebox.all()
        for i in shoebox:
            if i < 0:
                raise CentroidException, 'negative pixels in cube'
        data3d = shoebox.as_numpy_array()

        tot_itst = 0.0
        tot_f = 0.0
        tot_r = 0.0
        tot_c = 0.0

        for f in range(f_size):
            #print '__________________________________________________________________________________ new image'
            if numpy.sum(data3d[f, :, :]) > 0:
                locl_itst = single_spot_integrate_2d(data3d[f, :, :])
                tot_itst += locl_itst

        for f in range(f_size):
            for r in range(r_size):
                for c in range(c_size):
                    tot_f += f * data3d[f, r, c]
                    tot_r += r * data3d[f, r, c]
                    tot_c += c * data3d[f, r, c]

        if tot_itst > 0.0:
            _f = tot_f / tot_itst
            _r = tot_r / tot_itst
            _c = tot_c / tot_itst
        else:
            _f = -1
            _r = -1
            _c = -1
            raise CentroidException, 'negative or zero anything'

        #cont = 0
        tot_sf = 0.0
        tot_sr = 0.0
        tot_sc = 0.0

        for f in range(f_size):
            for r in range(r_size):
                for c in range(c_size):
                    tot_sf += ((f - _f) ** 2.0 * data3d[f, r, c])
                    tot_sr += ((r - _r) ** 2.0 * data3d[f, r, c])
                    tot_sc += ((c - _c) ** 2.0 * data3d[f, r, c])

        if tot_itst > 0.0:
            #_sf = numpy.sqrt(tot_sf) / tot_itst       # gaussian sigma formula
            #_sr = numpy.sqrt(tot_sr) / tot_itst       # gaussian sigma formula
            #_sc = numpy.sqrt(tot_sc) / tot_itst       # gaussian sigma formula

            _sf = tot_sf / tot_itst    # variance formula
            _sr = tot_sr / tot_itst    # variance formula
            _sc = tot_sc / tot_itst    # variance formula
        else:
            _sf = -1
            _sr = -1
            _sc = -1
        _f += 0.5
        _r += 0.5
        _c += 0.5

        return _f, _r, _c, _sf, _sr, _sc, tot_itst

def single_spot_integrate_2d(data2d):

    n_col = numpy.size(data2d[0:1, :])
    n_row = numpy.size(data2d[:, 0:1])
    #print 'n_col,n_row =', n_col, n_row

    data2dsmoth = numpy.zeros(n_row * n_col, dtype = float).reshape(n_row, n_col)

    diffdata2d = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    diffdata2d_ext = numpy.zeros(n_row * n_col, dtype = int).reshape(n_row, n_col)
    data2dtmp = data2d

    for times in range(5):
        for row in range(1, n_row - 1, 1):
            for col in range(1, n_col - 1, 1):
                pscan = float(numpy.sum(data2dtmp[row - 1:row + 2, col - 1:col + 2]) / 9.0)
                data2dsmoth[row, col] = pscan
        data2dtmp = data2dsmoth
    #threshold_shift = (numpy.max(data2dsmoth) - numpy.min(data2dsmoth)) / 2.0 # This used to be one of this "magical variables"

#######################################################################################################
    cont = 0                                                                  # This way to calculate
    dif_tot = 0                                                               # this magical variable
    for row in range(0, n_row, 1):                                            # is more statistical
        for col in range(0, n_col, 1):                                        # and seems to be giving
            cont += 1                                                         # better results
            dif_tot += numpy.abs(data2d[row, col] - data2dsmoth[row, col])    #
    dif_avg = dif_tot / cont                                                  #
    #print 'dif_avg=', dif_avg                                                #
    threshold_shift = dif_avg * 2.0                                           #
#######################################################################################################
    #print 'threshold_shift =', threshold_shift

    data2dsmoth[0:n_row, 0:n_col] = data2dsmoth[0:n_row, 0:n_col] + threshold_shift
    ext_area = 1                                                               # This used to be one of this "magical variables"
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            if data2d[row, col] > data2dsmoth[row, col]:
                diffdata2d[row, col] = 1

    for row in range(ext_area, n_row - ext_area + 1, 1):
        for col in range(ext_area, n_col - ext_area + 1, 1):
            if diffdata2d[row, col] == 1:
                diffdata2d_ext[row - ext_area:row + ext_area + 1, col - ext_area:col + ext_area + 1] = 1


###########################################################################flat background
    tot_bkgr = 0.0                                                       # version
    cont = 0.0                                                           #
    for row in range(0, n_row, 1):                                       #
        for col in range(0, n_col, 1):                                   #
            if diffdata2d_ext[row, col] == 0:                            #
                cont += 1                                                #
                tot_bkgr += data2d[row, col]                             #
    if tot_bkgr > 0 and cont > 0:                                        #
        bkgr = tot_bkgr / cont                                           #
    else:                                                                #
        bkgr = 0                                                         #
    #print 'bkgr=', bkgr                                                 #
    for row in range(0, n_row, 1):                                       #
        for col in range(0, n_col, 1):                                   #
            if diffdata2d_ext[row, col] == 1 and data2d[row, col] > bkgr:#
                data2d[row, col] = data2d[row, col] - bkgr               #
            else:                                                        #
                data2d[row, col] = 0                                     #
###########################################################################

######################################################################################## curved background
#   colbord = int(n_col / 5)                                                           # version
#   colbord = int(n_row / 5)                                                           #
#   for row in range(0, n_row, 1):                                                     #
#       for col in range(0, n_col, 1):                                                 #
#           if diffdata2d_ext[row,col] == 1:                                           #     #
#               top_av = float(numpy.sum(data2d[:colbord, col - 1:col + 2]))           #
#               bot_av = float(numpy.sum(data2d[n_row - colbord:, col - 1:col + 2]))   #
#               lft_av = float(numpy.sum(data2d[row - 1:row + 2, colbord]))            #
#               rgt_av = float(numpy.sum(data2d[row - 1:row + 2, n_col - colbord:]))   #
#               bkgr = (top_av + bot_av + lft_av + rgt_av) / 4.0                       #
#               if data2d[row,col] > bkgr:                                             #
#                   data2d[row,col] = data2d[row,col] - bkgr                           #
#               else:                                                                  #
#                   data2d[row,col] = 0                                                #
#   for row in range(0, n_row, 1):                                                     #
#       for col in range(0, n_col, 1):                                                 #
#           if diffdata2d_ext[row,col] == 0:                                           #
#               data2d[row,col] = 0                                                    #
###############################################################################

    itst_sum = 0.0
    for row in range(0, n_row, 1):
        for col in range(0, n_col, 1):
            itst_sum += float(data2d[row, col])
    if itst_sum > 0.0:
        tot = itst_sum
    else:
        print 'itst_sum =', itst_sum
        tot = -1

#    display_image_with_predicted_spots_n_centoids(data2d, col_cm, row_cm, n_col / 2, n_row / 2)

    return tot
def display_image_with_predicted_spots_n_centoids(image, xcoords, ycoords, xc, yc):
    """Display the image with coordinates overlayed."""
    from matplotlib import pylab, cm

    plt = pylab.imshow(image, vmin = 0, vmax = 1000, cmap = cm.Greys_r,
                       interpolation = 'nearest', origin = 'lower')
    pylab.scatter(xcoords, ycoords, marker = 'x')
    pylab.scatter(xc, yc, marker = 'x')
    plt.axes.get_xaxis().set_ticks([])
    plt.axes.get_yaxis().set_ticks([])
    pylab.show()
