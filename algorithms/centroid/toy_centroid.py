# toy centroid implementation: this is deliberately bad to illustrate how the
# interface works and how things could be improved
from __future__ import division
from dials.interfaces.centroid.centroid_interface_prototype import \
    centroid_interface_prototype as centroid_interface
from dials.interfaces.centroid.centroid_interface_prototype import \
    CentroidException

class toy_centroid(centroid_interface):
    def __init__(self, reflections):

        centroid_interface.__init__(self, reflections)

        return

    def compute_shoebox_centroid(self, shoebox):

        import math

        f_size, r_size, c_size = shoebox.all()

        # build the list of pixels - let's be dumb and just have a literal
        # list - and assign density of a pixel to the centre of the
        # volume...

        pixel_list = []

        for i in shoebox:
            if i < 0:
                raise CentroidException, 'negative pixels in cube'

        try:
            for f in range(f_size):
                for r in range(r_size):
                    for c in range(c_size):
                        pixel_list.append(
                            (f + 0.5, r + 0.5, c + 0.5, shoebox[f, r, c]))

        except IndexError, e:
            raise CentroidException, 'data outside range'

        # compute averages of positions

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0
        d_tot = 0.0

        for f, r, c, d in pixel_list:
            f_tot += d * f
            r_tot += d * r
            c_tot += d * c
            d_tot += d

        # print shoebox.as_numpy_array()

        if not d_tot:
            raise CentroidException("Invalid value for total intensity")

        _f, _r, _c = f_tot / d_tot, r_tot / d_tot, c_tot / d_tot

        # now compute the variance

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0

        for f, r, c, d in pixel_list:
            f_tot += d * (f - _f) ** 2
            r_tot += d * (r - _r) ** 2
            c_tot += d * (c - _c) ** 2

        _sf = f_tot / d_tot
        _sr = r_tot / d_tot
        _sc = c_tot / d_tot

        return _f, _r, _c, _sf, _sr, _sc, d_tot
