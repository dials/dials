
from dials.model.serialize import load
from dials.algorithms.reflection_basis import CoordinateSystem
from dials.algorithms.reflection_basis import CoordinateSystem
from dials.algorithms.reflection_basis import FromRotationAngleAccurate
from dials.algorithms.reflection_basis import ToRotationAngleAccurate
from dials.algorithms.reflection_basis import FromRotationAngleFast
from scitbx.array_family import flex
from scitbx import matrix
from math import pi, sqrt, atan
sweep = load.sweep('/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/sweep.json')

beam = sweep.get_beam()
gonio = sweep.get_goniometer()
detector = sweep.get_detector()
scan = sweep.get_scan()

m2 = matrix.col(gonio.get_rotation_axis())
s0 = matrix.col(beam.get_s0())
xsize, ysize = detector.get_image_size()
div = 4
mask1 = flex.bool(flex.grid(int(ysize / div), int(xsize / div)))
mask2 = flex.bool(flex.grid(int(ysize / div), int(xsize / div)))

n = 5
dm = -n * 0.157 * pi / 180.0
#dm = -n * 1.0 * pi / 180.0

#dm /= pi / 2

print dm

xyz = detector.get_pixel_lab_coord((0, 1200))
s1 = matrix.col(xyz).normalize() * s0.length()
phi = 0.0
cs = CoordinateSystem(m2, s0, s1, phi)
print "Limits: ", cs.limits()

#from_angle = FromRotationAngleAccurate(cs)
#phi_list = []
#e3_list = []
#for t in range(-100, 100):
#    p = phi - dm * t / 10
#    print p, from_angle(p)
#    phi_list.append(p)
#    e3_list.append(from_angle(p))

#from matplotlib import pylab
#pylab.plot(phi_list, e3_list)
#pylab.axvline(x=phi)
#pylab.show()



m2 = matrix.col(cs.m2())
e1 = matrix.col(cs.e1_axis())
e3 = matrix.col(cs.e3_axis())
ps = matrix.col(cs.p_star()).normalize()

m2e1 = m2.dot(e1)
m2e3 = m2.dot(e3)
m2ps = m2.dot(ps)
ev = m2e1**2 + 2*dm*m2e3*m2ps - dm**2

inc1 = ev >= 0.0
zeta = cs.zeta()

inc2 = abs(zeta) >= 0.05
print ev, zeta

a = m2e1
b = m2e3*m2ps
tanphi0 = (b + sqrt(a**2 + b**2)) / a
tanphi1 = (b - sqrt(a**2 + b**2)) / a
print 2.0 * atan(tanphi0), 2.0 * atan(tanphi1)

print 1/0





for j in range(int(ysize / div)):
    print j
    for i in range(int(xsize / div)):

        xyz = detector.get_pixel_lab_coord((div*i, div*j))
        s1 = matrix.col(xyz).normalize() * s0.length()
        phi = 0.0
        cs = CoordinateSystem(m2, s0, s1, phi)
        m2 = matrix.col(cs.m2())
        e1 = matrix.col(cs.e1_axis())
        e3 = matrix.col(cs.e3_axis())
        ps = matrix.col(cs.p_star()).normalize()

        from_angle = FromRotationAngle(cs)
        c3 = from_angle(dm)

        m2e1 = m2.dot(e1)
        m2e3 = m2.dot(e3)
        m2ps = m2.dot(ps)
        ev = m2e1**2 + 2*dm*m2e3*m2ps - dm**2

        inc1 = ev >= 0.0
        zeta = cs.zeta()

#        if div*i == 0 and div*j == 600:
#            print ev, zeta
#            print 1/0

        inc2 = abs(zeta) >= 0.05

        mask1[j,i] = inc1
        mask2[j,i] = inc2

from matplotlib import pylab, cm
pylab.subplot(1, 2, 1)
pylab.imshow(mask1.as_numpy_array(), cmap=cm.Greys_r)
pylab.subplot(1, 2, 2)
pylab.imshow(mask2.as_numpy_array(), cmap=cm.Greys_r)
pylab.show()
