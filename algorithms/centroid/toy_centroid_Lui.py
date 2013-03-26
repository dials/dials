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

#        import math

        f_size, r_size, c_size = shoebox.all()


        for i in shoebox:
            if i < 0:
                raise CentroidException, 'negative pixels in cube'
        data3d = shoebox.as_numpy_array()

        tot_itst = 0.0
        tot_f = 0.0
        tot_c = 0.0
        tot_r = 0.0
        tot_sr = 0.0
        tot_sc = 0.0

        cont = 0
        itst = numpy.zeros(f_size, dtype = float).reshape(f_size)
        #for f in range(f_min, f_max):
        for f in range(f_size):
            print '__________________________________________________________________________________ new image'
            data2d = data3d[f, :, :]
            if numpy.sum(data2d) > 0:
                print '______________________data2d before'
                print data2d
                row_cm, col_cm, locl_sr, locl_sc, locl_itst = single_spot_integrate_2d(data2d)
                print '______________________data2d after'
                print data2d

                itst[cont] = locl_itst
                cont += 1
                tot_f += f * locl_itst
                tot_r += row_cm * locl_itst
                tot_c += col_cm * locl_itst
                tot_itst += locl_itst
                #print f, locl_itst, locl_f
                tot_sr += (locl_sr * locl_itst) ** 2.0
                tot_sc += (locl_sc * locl_itst) ** 2.0

        if tot_itst > 0.0:
            _f = tot_f / tot_itst
            _r = tot_r / tot_itst
            _c = tot_c / tot_itst
        else:
            _f = -1
            _r = -1
            _c = -1
        cont = 0
        tot_sf = 0.0

        #for f in range(f_min, f_max):
        for f in range(f_size):
            tot_sf += ((f - _f) * itst[cont]) ** 2.0
            cont += 1
        if tot_itst > 0.0:
            _sf = numpy.sqrt(tot_sf) / tot_itst
            _sr = numpy.sqrt(tot_sr) / tot_itst
            _sc = numpy.sqrt(tot_sc) / tot_itst
        else:
            _sf = -1
            _sr = -1
            _sc = -1
        print '_f, _r, _c, _sf, _sr, _sc =', _f, _r, _c, _sf, _sr, _sc

        return _f, _r, _c, _sf, _sr, _sc

def single_spot_integrate_2d(data2d):

    x_to = numpy.size(data2d[0:1, :])
    y_to = numpy.size(data2d[:, 0:1])
    print 'x_to,y_to =', x_to, y_to

    data2dsmoth = numpy.zeros(y_to * x_to, dtype = float).reshape(y_to, x_to)

    diffdata2d = numpy.zeros(y_to * x_to, dtype = int).reshape(y_to, x_to)
    diffdata2d_ext = numpy.zeros(y_to * x_to, dtype = int).reshape(y_to, x_to)
    data2dtmp = data2d

    ext_area = 1                                                               # This used to be one of this "magical variables"

    for times in range(5):
        for y in range(1, y_to - 1, 1):
            for x in range(1, x_to - 1, 1):
                pscan = float(numpy.sum(data2dtmp[y - 1:y + 2, x - 1:x + 2]) / 9.0)
                data2dsmoth[y, x] = pscan
    data2dtmp = data2dsmoth
    #threshold_shift = (numpy.max(data2dsmoth) - numpy.min(data2dsmoth)) / 2.0 # This used to be one of this "magical variables"

#######################################################################################################
    cont = 0                                                           # This way to calculate
    dif_tot = 0                                                        # this magical variable
    for y in range(0, y_to, 1):                                        # is more statistical
        for x in range(0, x_to, 1):                                    # and seems to be giving
            cont += 1                                                  # better results
            dif_tot += numpy.abs(data2d[y, x] - data2dsmoth[y, x])     #
    dif_avg = dif_tot / cont                                           #
    print 'dif_avg=', dif_avg                                          #
    threshold_shift = dif_avg * 2.0                                    #
#######################################################################################################
    print 'threshold_shift =', threshold_shift

    data2dsmoth[0:y_to, 0:x_to] = data2dsmoth[0:y_to, 0:x_to] + threshold_shift

    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if data2d[y, x] > data2dsmoth[y, x]:
                diffdata2d[y, x] = 1

    for y in range(ext_area, y_to - ext_area + 1, 1):
        for x in range(ext_area, x_to - ext_area + 1, 1):
            if diffdata2d[y, x] == 1:
                diffdata2d_ext[y - ext_area:y + ext_area + 1, x - ext_area:x + ext_area + 1] = 1

    print '_____________diffdata2d'
    #print diffdata2d
    print '_____________diffdata2d_ext'
    #print diffdata2d_ext


############################################################################### flat background
    tot_bkgr = 0.0                                                            # version
    cont = 0.0                                                                #
    for y in range(0, y_to, 1):                                               #
        for x in range(0, x_to, 1):                                           #
            if diffdata2d_ext[y, x] == 0:                                     #
                cont += 1                                                     #
                tot_bkgr += data2d[y, x]                                      #
    bkgr = tot_bkgr / cont                                                    #
    print 'bkgr=', bkgr                                                       #
    for y in range(0, y_to, 1):                                               #
        for x in range(0, x_to, 1):                                           #
            if diffdata2d_ext[y, x] == 1 and data2d[y, x] > bkgr:             #
                data2d[y, x] = data2d[y, x] - bkgr                            #
            else:                                                             #
                data2d[y, x] = 0                                              #
###############################################################################

############################################################################### curved background
#   xbord = int(x_to / 5)                                                     # version
#   ybord = int(y_to / 5)                                                     #
#   for y in range(0, y_to, 1):                                               #
#       for x in range(0, x_to, 1):                                           #
#           if diffdata2d_ext[y, x] == 1:                                     #
#               top_av = float(numpy.sum(data2d[:ybord, x - 1:x + 2]))        #
#               bot_av = float(numpy.sum(data2d[y_to - ybord:, x - 1:x + 2])) #
#               lft_av = float(numpy.sum(data2d[y - 1:y + 2, xbord]))         #
#               rgt_av = float(numpy.sum(data2d[y - 1:y + 2, x_to - xbord:])) #
#               bkgr = (top_av + bot_av + lft_av + rgt_av) / 4.0              #
#               if data2d[y, x] > bkgr:                                       #
#                   data2d[y, x] = data2d[y, x] - bkgr                        #
#               else:                                                         #
#                   data2d[y, x] = 0                                          #
#   for y in range(0, y_to, 1):                                               #
#       for x in range(0, x_to, 1):                                           #
#           if diffdata2d_ext[y, x] == 0:                                     #
#               data2d[y, x] = 0                                              #
###############################################################################
    x_num_sum = 0.0
    y_num_sum = 0.0
    den_sum = 0.0
    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            x_num_sum = x_num_sum + float(data2d[y, x]) * float(x)
            y_num_sum = y_num_sum + float(data2d[y, x]) * float(y)
            den_sum = den_sum + float(data2d[y, x])
    if den_sum > 0.0:
        col_cm = x_num_sum / den_sum
        row_cm = y_num_sum / den_sum
    else:
        print 'den_sum =', den_sum
        col_cm = -1
        row_cm = -1

    x_num_sum = 0.0
    y_num_sum = 0.0
#    den_sum = 0.0
    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            x_num_sum = x_num_sum + float(data2d[y, x]) * (float(x) - float(col_cm)) ** 2.0
            y_num_sum = y_num_sum + float(data2d[y, x]) * (float(y) - float(row_cm)) ** 2.0
    #        den_sum = den_sum + float(data2d[y, x])
    if den_sum > 0.0:
        col_sig = numpy.sqrt(x_num_sum / den_sum)
        row_sig = numpy.sqrt(y_num_sum / den_sum)
    else:
        print 'den_sum =', den_sum
        col_sig = -1
        row_sig = -1

#    display_image_with_predicted_spots_n_centoids(data2d, col_cm, row_cm, x_to / 2, y_to / 2)

    return row_cm, col_cm, row_sig, col_sig, den_sum

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
