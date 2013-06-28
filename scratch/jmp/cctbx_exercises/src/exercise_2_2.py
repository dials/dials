from __future__ import division

from scitbx import matrix

h = -54
k = 9
l = 5
s0 = matrix.col((0.013141, 0.002201, 1.450476))
m2 = matrix.col((0.999975, -0.001289, -0.006966))
Astar = matrix.sqr((33.914,-36.200, -14.128,
                38.856, 31.328, 13.000,
                -0.543, -19.192, 47.871)).inverse()
z0 = 1
z = 177.1
phi0 = 1.0000
dphi = 1.0000
psi = -144.87

phi = phi0 + (z - z0) * dphi

print "Phi: ", phi

p_star0 = Astar * (h, k, l)

p_star = p_star0.rotate(axis=m2, angle=phi, deg=True)
s = s0 + p_star

print "S0 length:", s0.length()
print "S length:", s.length()

print 1/0



from rstbx.diffraction import rotation_angles
from scitbx import matrix
import math

hkl = matrix.col((h, k, l))
wavelength = 0.689400
dmin = 0.7 # guessed
ub_matrix = matrix.sqr(tuple(b1) + tuple(b2) + tuple(b3)).inverse()

ra = rotation_angles(dmin, ub_matrix, wavelength, matrix.col((1, 0, 0)))
r2d = 180.0 / math.pi
if ra(hkl):
    phis = ra.get_intersection_angles()
    phical0 = phis[0] * r2d % 360
    phical1 = phis[1] * r2d % 360
    print "Predicted Phi:", (phical0, hkl)
    print "Predicted Phi:", (phical1, hkl)
    zcal0 = (phical0 - phi0) / dphi + z0
    zcal1 = (phical1 - phi0) / dphi + z0
    print "Predicted Frames:", zcal0, zcal1

    print phical0 - phi
    print phical1 - phi
