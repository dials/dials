# toy centroid implementation: this is deliberately bad to illustrate how the
# interface works and how things could be improved
from __future__ import division
from dials.interfaces.centroid.centroid_interface_prototype import \
    centroid_interface_prototype as centroid_interface

def covar_contrib(x1, x2, weights):
    '''Compute contribution to covariance matrix.'''

    covar = sum([_x1 * _x2 * _w for _x1, _x2, _w in zip(x1, x2, weights)]) / \
        sum(weights)

    return covar
        

class toy_centroid(centroid_interface):
    def __init__(self, bounding_boxes, dxtbx_sweep_object):

        self._image_size = dxtbx_sweep_object.get_detector().get_image_size()

        centroid_interface.__init__(self, bounding_boxes, dxtbx_sweep_object)


        return

    def compute_centroid_from_bbox(self, bbox):

        import math

        f_min, f_max, r_min, r_max, c_min, c_max = bbox

        # build the list of pixels - let's be dumb and just have a literal
        # list

        pixel_list = []

        try:
            for f in range(f_min, f_max):
                data = self._sweep[f]
                for r in range(r_min, r_max):
                    for c in range(c_min, c_max):
                        pixel_list.append(
                            (f, r, c, data[r * self._image_size[0] + c]))

        except IndexError, e:
            return -1., -1., -1., -1., -1., -1.

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

        assert(d_tot)

        _f, _r, _c = f_tot / d_tot, r_tot / d_tot, c_tot / d_tot

        # form the variance matrix

        pixels = [(f - _f, r - _r, c - _c, d) for f, r, c, d in pixel_list]

        data = ([pixel[0] for pixel in pixels],
                [pixel[1] for pixel in pixels],
                [pixel[2] for pixel in pixels],
                [pixel[3] for pixel in pixels])

        m_elems = []

        for i in range(3):
            for j in range(3):
                m_elems.append(covar_contrib(data[i], data[j], data[3]))

        from scitbx.array_family import flex
        from scitbx.linalg import eigensystem
        from scitbx import matrix

        m = flex.double(flex.grid(3, 3))
        for j in range(9):
            m[j] = m_elems[j]
        s = eigensystem.real_symmetric(m)
        values = s.values()
        vectors = s.vectors()

        print 'Eigenvalues and vectors of covariance matrix from spot...'
        for j in range(3):
            vec = tuple(vectors[j * 3:j * 3 + 3])
            print '%.3f %.3f %.3f %.3f' % (values[j], vec[0], vec[1], vec[2])
            
        # now compute the variance

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0
        d_tot = 0.0

        for f, r, c, d in pixel_list:
            f_tot += d * (f - _f) ** 2
            r_tot += d * (r - _r) ** 2
            c_tot += d * (c - _c) ** 2
            d_tot += d

        _sf = math.sqrt(f_tot / d_tot)
        _sr = math.sqrt(r_tot / d_tot)
        _sc = math.sqrt(c_tot / d_tot)

        return _f, _r, _c, _sf, _sr, _sc
