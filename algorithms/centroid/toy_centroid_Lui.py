# toy centroid implementation: this is deliberately bad to illustrate how the
# interface works and how things could be improved
# from __future__ import division
from dials.interfaces.centroid.centroid_interface_prototype import \
    centroid_interface_prototype as centroid_interface
import numpy

class toy_centroid_lui(centroid_interface):
    def __init__(self, bounding_boxes, dxtbx_sweep_object):

        self._image_size = dxtbx_sweep_object.get_detector().get_image_size()

        centroid_interface.__init__(self, bounding_boxes, dxtbx_sweep_object)

        return

    def compute_centroid_from_bbox(self, bbox):

        import math

        f_min, f_max, r_min, r_max, c_min, c_max = bbox

        # build a numpy array that represents a 3D block of pixels
        f_size = f_max - f_min
        r_size = r_max - r_min
        c_size = c_max - c_min
        data3d = numpy.zeros(f_size * r_size * c_size, dtype = int).reshape(f_size, r_size, c_size)
        in_f = 0
        in_r = 0
        in_c = 0
        for f in range(f_min, f_max):
            data = self._sweep[f]
            for r in range(r_min, r_max):
                for c in range(c_min, c_max):
                    data3d[in_f, in_r, in_c] = data[r * self._image_size[0] + c]
                    in_c += 1
                in_c = 0
                in_r += 1
            in_r = 0
            in_f += 1

        tot_itst = 0.0
        tot_f = 0.0
        tot_c = 0.0
        tot_r = 0.0
        tot_sr = 0.0
        tot_sc = 0.0

        cont = 0
        itst = numpy.zeros(f_size, dtype = float).reshape(f_size)
        for f in range(f_min, f_max):
            data2d = data3d[f, :, :]
            print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            print data2d
            row_cm, col_cm, locl_sr, locl_sc, locl_itst = single_spot_integrate_2d(data2d)
            print data2d
            print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            itst[cont] = locl_itst
            cont += 1
            locl_f = f * locl_itst
            tot_f += locl_f
            tot_r += (row_cm + r_min) * locl_itst
            tot_c += (col_cm + c_min) * locl_itst
            tot_itst += locl_itst
            #print f, locl_itst, locl_f
            tot_sr += (locl_sr * locl_itst) ** 2.0
            tot_sc += (locl_sc * locl_itst) ** 2.0

        _f = tot_f / tot_itst
        _r = tot_r / tot_itst
        _c = tot_c / tot_itst
        cont = 0
        tot_sf = 0.0

        for f in range(f_min, f_max):
            tot_sf += ((f - _f) * itst[cont]) ** 2.0
            cont += 1

        _sf = numpy.sqrt(tot_sf) / tot_itst
        _sr = numpy.sqrt(tot_sr) / tot_itst
        _sc = numpy.sqrt(tot_sc) / tot_itst

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

    ext_area = 2                                                               # This used to be one of this "magical variables"

    for times in range(5):
        for y in range(1, y_to - 1, 1):
            for x in range(1, x_to - 1, 1):
                pscan = float(numpy.sum(data2dtmp[y - 1:y + 2, x - 1:x + 2]) / 9.0)
                data2dsmoth[y, x] = pscan
    data2dtmp = data2dsmoth
    threshold_shift = numpy.ptp(data2dsmoth)                                   # This used to be one of this "magical variables"
    print 'threshold_shift =', threshold_shift
    data2dsmoth[0:y_to, 0:x_to] = data2dsmoth[0:y_to, 0:x_to] + threshold_shift

    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if data2d[y, x] > data2dsmoth[y, x]:
                diffdata2d[y, x] = 1

    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if diffdata2d[y, x] == 1:
                diffdata2d_ext[y - ext_area:y + ext_area + 1, x - ext_area:x + ext_area + 1] = 1
    xbord = int(x_to / 5)
    ybord = int(y_to / 5)
    print 'xbord, ybord =', xbord, ybord

    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if diffdata2d_ext[y, x] == 1:
                top_av = float(numpy.sum(data2d[:ybord, x - 1:x + 2]))
                bot_av = float(numpy.sum(data2d[y_to - ybord:, x - 1:x + 2]))
                lft_av = float(numpy.sum(data2d[y - 1:y + 2, xbord]))
                rgt_av = float(numpy.sum(data2d[y - 1:y + 2, x_to - xbord:]))
                bkgr = (top_av + bot_av + lft_av + rgt_av) / 4.0
                if data2d[y, x] > bkgr:
                    data2d[y, x] = data2d[y, x] - bkgr
    for y in range(0, y_to, 1):
        for x in range(0, x_to, 1):
            if diffdata2d_ext[y, x] == 0:
                data2d[y, x] = 0

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
