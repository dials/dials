def lorentz(m2, s0, s1):
    return abs(m2.dot(s1.cross(s0))) / (s1.length() * s0.length())

def include_ref(m2, e1, e3, ps, dm):
    m2e1 = m2.dot(e1)
    m2e3 = m2.dot(e3)
    m2ps = m2.dot(ps)
    return dm**2 - 2 * dm * m2e3 * m2ps - m2e3**2 >= 0.0


from math import pi, acos
from scitbx import matrix
from dials.algorithms.reflection_basis import CoordinateSystem
from dials.algorithms.reflection_basis import FromRotationAngleAccurate
from dials.algorithms.reflection_basis import ToRotationAngleAccurate
from dials.algorithms.reflection_basis import FromRotationAngleFast

#s0 = matrix.col((0.00801639379479156, -3.8514801707506e-14, -1.0208975723528189) )
#m2 = matrix.col((1.0, -1.5919306617286774e-16, -6.904199434387693e-16))
#s1 = matrix.col((0.502248600762, -0.0181543707198, -0.888657908118))

#s0 = matrix.col((0.0, 0.0, -1.0))
#m2 = matrix.col((1.0, 0.0, 0.0))
#s1 = matrix.col((1.0, 0.05, 0.0)).normalize()
phi =  0

#c3 = -0.013439
n = 5
c3 = -n * 0.157 * pi / 180.0
cs = CoordinateSystem(m2, s0, s1, phi)
print cs.zeta()
#e1 = matrix.col(cs.e1_axis())
#e3 = matrix.col(cs.e3_axis())
#ps = matrix.col(cs.p_star())

#from math import sqrt
#m2e1 = m2.dot(e1)
#m2e3 = m2.dot(e3)
#m2ps = m2.dot(ps.normalize())

#print tuple(e1)
#print tuple(e3)
#print tuple(ps)
#print "m2.e1", m2e1
#print "m2.e3", m2e3
#print "m2.ps", m2ps
#print "(m2.e3)(m2.ps)", m2e3*m2ps, 1.0 + m2e1
#print "Angle: ", s0.angle(s1, deg=True)
#print "Limits: ", cs.limits()
#print "Lorentz: ", lorentz(m2, s0, s1)
#print m2.length(), e1.length(), e3.length()
#print "Length: ", sqrt(m2e3**2 + m2ps**2 + m2e1**2)
#print "Sum Angles: ", (acos(m2e3) + acos(m2ps) + acos(m2e1)) * 180 / pi
#to_angle = ToRotationAngleAccurate(cs)
#phi = to_angle(cs.limits()[2])



#print phi

from_angle = FromRotationAngleAccurate(cs)
phi_list = []
e3_list = []
for t in range(-10, 10):
    p = phi - 0.1 * t * pi / 180
    print p, from_angle(p)
    phi_list.append(p * 180 / pi)
    e3_list.append(from_angle(p))

from matplotlib import pylab
pylab.plot(phi_list, e3_list)
pylab.axvline(x=phi)
pylab.show()

#to_phi = rb.ToRotationAngleAccurate(cs)
#m2e1 = matrix.col(cs.m2()).dot(matrix.col(cs.e1_axis()))
#m2e1_m2e1 = m2e1*m2e1
#m2e3_m2ps = 2.0 * (matrix.col(cs.m2()).dot(matrix.col(cs.e3_axis()))) *\
#                  (matrix.col(cs.m2()).dot(matrix.col(cs.p_star()).normalize()))
#print m2e1, m2e1_m2e1, m2e3_m2ps
#print m2e1_m2e1 + c3 * m2e3_m2ps - c3*c3
#print to_phi(-0.013439), to_phi(0.013439)
