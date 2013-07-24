
from dials.model.serialize import load
from dials.algorithms.reflection_basis import CoordinateSystem
from scitbx.array_family import flex
from scitbx import matrix
from math import pi
sweep = load.sweep('/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/sweep.json')

beam = sweep.get_beam()
gonio = sweep.get_goniometer()
detector = sweep.get_detector()
scan = sweep.get_scan()

m2 = matrix.col(gonio.get_rotation_axis())
s0 = matrix.col(beam.get_s0())
xsize, ysize = detector.get_image_size()

div = 4
mask1 = flex.bool(flex.grid(int(ysize / 4), int(xsize / 4)))
mask2 = flex.bool(flex.grid(int(ysize / 4), int(xsize / 4)))

n = 5
#dm = -n * 0.157 * pi / 180.0
dm = -n * 1.0 * pi / 180.0

for j in range(int(ysize / 4)):
    print j
    for i in range(int(xsize / 4)):

        xyz = detector.get_pixel_lab_coord((4*i, 4*j))
        s1 = matrix.col(xyz).normalize() * s0.length()
        phi = 0.0
        cs = CoordinateSystem(m2, s0, s1, phi)
        e1 = matrix.col(cs.e1_axis())
        e3 = matrix.col(cs.e3_axis())
        ps = matrix.col(cs.p_star()).normalize()

        m2e1 = m2.dot(e1)
        m2e3 = m2.dot(e3)
        m2ps = m2.dot(ps)
        ev = m2e1**2 + 2*dm*m2e3*m2ps - dm**2

        inc1 = ev >= 0.0
        zeta = cs.zeta()

        inc2 = abs(zeta) >= 0.5

        mask1[j,i] = inc1
        mask2[j,i] = inc2

from matplotlib import pylab, cm
pylab.subplot(1, 2, 1)
pylab.imshow(mask1.as_numpy_array(), cmap=cm.Greys_r)
pylab.subplot(1, 2, 2)
pylab.imshow(mask2.as_numpy_array(), cmap=cm.Greys_r)
pylab.show()
