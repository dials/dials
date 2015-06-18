#!/usr/bin/env dials.python
#
# Play area trying to derive useful derivatives for Rodrigues rotation formula:
#
# http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
#
# Details including sage script in dx.tex

def skew_symm(v):
  '''Make matrix [v]_x from v. Essentially multiply vector by SO(3) basis
  set Lx, Ly, Lz. Equation (2) from Gallego & Yezzi paper.'''
  import scitbx.matrix

  L1 = scitbx.matrix.sqr((0, 0, 0, 0, 0, -1, 0, 1, 0))
  L2 = scitbx.matrix.sqr((0, 0, 1, 0, 0, 0, -1, 0, 0))
  L3 = scitbx.matrix.sqr((0, -1, 0, 1, 0, 0, 0, 0, 0))

  v1, v2, v3 = v.elems

  return v1 * L1 + v2 * L2 + v3 * L3

def dR_drs_fd(r, s, t):
  # perform finite difference calculation of dR/dk1 as
  # (1/2 eps) * (R(k1+eps) - R(k1-eps))

  eps = 1.0e-8

  Rpls = Rtk(t, krs(r + eps, s))
  Rmns = Rtk(t, krs(r - eps, s))
  ndRr = (0.5 / eps) * (Rpls - Rmns)

  Rpls = Rtk(t, krs(r, s + eps))
  Rmns = Rtk(t, krs(r, s - eps))
  ndRs = (0.5 / eps) * (Rpls - Rmns)

  return ndRr, ndRs

def dR_dki(t, k):
  '''This is equation (7) from Gallego & Yezzi'''
  from scitbx import matrix
  from math import sin, cos

  K = skew_symm(k)

  e = (matrix.col((1, 0, 0)),
       matrix.col((0, 1, 0)),
       matrix.col((0, 0, 1)))

  return [
    t * (cos(t) * k[i] * K + sin(t) * k[i] * K * K +
      (sin(t)/t) * skew_symm(e[i] - k[i] * k) +
      ((1 - cos(t)) / t) * (e[i] * k.transpose() + k * e[i].transpose() -
                            2 * k[i] * k * k.transpose())) for i in range(3)]

def gallego_yezzi_eqn9(t, k):
  '''This is equation (9) from Gallego & Yezzi'''

  from scitbx import matrix
  from math import sin, cos

  v = t*k
  R = k.axis_and_angle_as_r3_rotation_matrix(t, deg=False)
  I3 = matrix.identity(3)

  V = skew_symm(v)

  e = (matrix.col((1, 0, 0)),
       matrix.col((0, 1, 0)),
       matrix.col((0, 0, 1)))

  term1 = [skew_symm(v.cross((I3 - R) * e[i])) for i in range(3)]

  return [
    (1./t) * (v[i] * V + term1[i]) * R for i in range(3)]


def dR_drs_calc(r, s, t):
  from scitbx import matrix
  from math import sin, cos

  # pre-compile some useful values (makes code which follows slightly easier
  # to read

  cr = cos(r); sr = sin(r); cs = cos(s); ss = sin(s); ct = cos(t); st = sin(t)
  cr2 = cos(r * 2); sr2 = sin(r * 2); cs2 = cos(s * 2); ss2 = sin(s * 2)

  dRdr = matrix.sqr([
    -1/2*(cs2*sr2 - sr2)*ct + 1/2*cs2*sr2 - cr*sr,
    cr**2 - 1/2*cr2*cs2 + 1/2*(cr2*cs2 - cr2)*ct - 1/2,
    1/2*ct*sr*ss2 + cr*ss*st - 1/2*sr*ss2,
    cr**2 - 1/2*cr2*cs2 + 1/2*(cr2*cs2 - cr2)*ct - 1/2,
    1/2*(cs2*sr2 - sr2)*ct - 1/2*cs2*sr2 + cr*sr,
    -1/2*cr*ct*ss2 + sr*ss*st + 1/2*cr*ss2,
    1/2*ct*sr*ss2 - cr*ss*st - 1/2*sr*ss2,
    -1/2*cr*ct*ss2 - sr*ss*st + 1/2*cr*ss2,
    0])

  dRds = matrix.sqr([
    -1/2*(cr2 + 1)*ct*ss2 + 1/2*cr2*ss2 + cs*ss,
    -1/2*ct*sr2*ss2 + 1/2*sr2*ss2 + ss*st,
    -cr*cs2*ct + cs*sr*st + cr*cs2,
    -1/2*ct*sr2*ss2 + 1/2*sr2*ss2 - ss*st,
    1/2*(cr2 - 1)*ct*ss2 - 1/2*cr2*ss2 + cs*ss,
    -cs2*ct*sr - cr*cs*st + cs2*sr,
    -cr*cs2*ct - cs*sr*st + cr*cs2,
    -cs2*ct*sr + cr*cs*st + cs2*sr,
    ct*ss2 - 2*cs*ss])

  return dRdr, dRds

def Rtk(t, k):
  import scitbx.matrix
  import math

  I = scitbx.matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  c = math.cos(t)
  s = math.sin(t)
  K = skew_symm(k)

  R = I + s * K + (1 - c) * K * K

  # make sure answer agrees with scitbx matrix calculations (optional)
  # zero = R.inverse() * k.axis_and_angle_as_r3_rotation_matrix(t) - I
  # assert sum(zero.elems) < 1.0e-8

  return R

def krs(r, s):
  from math import sin, cos
  import scitbx.matrix
  k1 = cos(r) * sin(s)
  k2 = sin(r) * sin(s)
  k3 = cos(s)
  return scitbx.matrix.col((k1, k2, k3))

def dk_dr(r, s):
  from math import sin, cos
  import scitbx.matrix
  k1 = -sin(r) * sin(s)
  k2 = cos(r) * sin(s)
  k3 = 0
  return scitbx.matrix.col((k1, k2, k3))

def dk_ds(r, s):
  from math import sin, cos
  import scitbx.matrix
  k1 = cos(r) * cos(s)
  k2 = sin(r) * cos(s)
  k3 = -sin(s)
  return scitbx.matrix.col((k1, k2, k3))

def work():
  import random
  import math

  # pick random rotation axis & angle - pick spherical coordinate basis to make
  # life easier in performing finite difference calculations

  t = 2 * math.pi * random.random()
  r = 2 * math.pi * random.random()
  s = math.pi * random.random()

  # finite difference version
  ndRr, ndRs = dR_drs_fd(r, s, t)

  # now call on derivative calculation
  dRr, dRs = dR_drs_calc(r, s, t)

  assert sum((ndRr - dRr).elems) < 1.0e-7
  assert sum((ndRs - dRs).elems) < 1.0e-7

def dx():
  for j in range(100):
    work()
  print 'OK'

# dx()

def work():
  import random
  import math
  from scitbx import matrix

  # pick random rotation axis & angle - pick spherical coordinate basis to make
  # life easier in performing finite difference calculations

  t = 2 * math.pi * random.random()
  r = 2 * math.pi * random.random()
  s = math.pi * random.random()
  k = krs(r, s)

  # finite difference version
  ndRr, ndRs = dR_drs_fd(r, s, t)

  # now call on derivative calculation
  dRr, dRs = dR_drs_calc(r, s, t)

  # G&Y eqn 7
  dRk = dR_dki(t, k)

  # G&Y eqn 9
  dRk2 = gallego_yezzi_eqn9(t, k)

  dkdr = dk_dr(r, s)
  dkds = dk_ds(r, s)

  print 'dR/dr'
  print 'FD'
  print ndRr
  print 'Analytical (G&Y 7)'
  m = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
  for i in range(3):
    m += dkdr[i] * dRk[i]
  print m
  print 'Analytical (G&Y 9)'
  m = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
  for i in range(3):
    m += dkdr[i] * dRk2[i]
  print m
  print 'Analytical (old)'
  print dRr

  print 'dR/ds'
  print 'FD'
  print ndRs
  print 'Analytical (G&Y 7)'
  m = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
  for i in range(3):
    m += dkds[i] * dRk[i]
  print m
  print 'Analytical (G&Y 9)'
  m = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
  for i in range(3):
    m += dkds[i] * dRk2[i]
  print m
  print 'Analytical (old)'
  print dRs

  # Finally test the direct derivative of a rotated vector wrt the axis elements
  # versus finite differences

  # make a random vector to rotate
  u = matrix.col((
    random.random(),
    random.random(),
    random.random()))

  # calc derivatives of rotated vector (Gallego & Yezzi equn 8)
  from dials_refinement_helpers_ext import dRq_de
  from scitbx.array_family import flex
  dr_de = dRq_de(flex.double([t]), flex.vec3_double([k]), flex.vec3_double([u]))
  print
  print "d[r]/d[e], where [r] = [R][u] is a rotation about [e] (G&Y 8)"
  print matrix.sqr(dr_de[0])

  print "Compare with FD calculation"
  dr_de_FD = [dR_ki * u for dR_ki in dRk]
  dr_de_FD = [elt for vec in dr_de_FD for elt in vec] # flatten list
  dr_de_FD = matrix.sqr(dr_de_FD).transpose() # make a matrix
  print dr_de_FD


work()

def tst_use_in_stills_parameterisation(beam_param=0):

  # test use of analytical expression in stills prediction parameterisation

  from scitbx import matrix
  from math import pi
  import random

  print
  print "Test use of analytical expressions in stills prediction parameterisation"

  # beam model
  from dxtbx.model.experiment import beam_factory
  from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisation
  beam = beam_factory().make_beam(matrix.col((0,0,1)), wavelength=1.1)
  s0 = matrix.col(beam.get_s0())
  s0u = matrix.col(beam.get_unit_s0())

  # beam parameterisation
  bp = BeamParameterisation(beam)

  # choose the second parameter, mu2, which is the one we usually refine. The
  # derivative of this is a vector along the X direction.
  ds0_dp = bp.get_ds_dp()[beam_param]

  # pick a random point on (the positive octant of) the Ewald sphere to rotate
  s1 = matrix.col((
    random.random(),
    random.random(),
    random.random())).normalize() * s0.length()
  r = s1 - s0
  r0 = r.normalize()

  # calculate the axis of rotation
  e1 = r0.cross(s0u).normalize()

  # calculate c0, a vector orthogonal to s0u and e1
  c0 = s0u.cross(e1)

  # convert to derivative of the unit beam direction. This involves scaling
  # by the wavelength, then projection back onto the Ewald sphere.
  scaled = ds0_dp * beam.get_wavelength()
  ds0u_dp = scaled.dot(c0) * c0 + scaled.dot(e1) * e1

  # rotate relp off Ewald sphere a small angle (up to 1 deg)
  DeltaPsi = random.uniform(-pi/180, pi/180)
  q = r.rotate_around_origin(e1, -DeltaPsi)
  q0 = q.normalize()
  from libtbx.test_utils import approx_equal
  assert approx_equal(q0.cross(s0u).normalize(), e1)

  # use the fact that e1 == c0.cross(s0u) to redefine the derivative d[e1]/dp
  # from Sauter et al. (2014) (A.3)
  de1_dp = c0.cross(ds0u_dp)

  # unlike the previous definition this *is* orthogonal to e1, as expected.
  print "[e1].(d[e1]/dp) = {0} (should be 0.0)".format(e1.dot(de1_dp))

  # calculate (d[r]/d[e1])(d[e1]/dp) analytically
  from scitbx.array_family import flex
  from dials_refinement_helpers_ext import dRq_de
  dr_de1 = matrix.sqr(dRq_de(flex.double([DeltaPsi]),
                  flex.vec3_double([e1]),
                  flex.vec3_double([q]))[0])
  print "Analytical calculation for (d[r]/d[e1])(d[e1]/dp):"
  print dr_de1 * de1_dp

  # now calculate using finite differences.
  dp = 1.e-8
  del_e1 = de1_dp * dp
  e1f = e1 + del_e1 * 0.5
  rfwd = q.rotate_around_origin(e1f, DeltaPsi)
  e1r = e1 - del_e1 * 0.5
  rrev = q.rotate_around_origin(e1r, DeltaPsi)

  print "Finite difference estimate for (d[r]/d[e1])(d[e1]/dp):"
  print (rfwd - rrev) * (1 / dp)

  print "These are now the same :-)"

tst_use_in_stills_parameterisation(0)
tst_use_in_stills_parameterisation(1)
tst_use_in_stills_parameterisation(2)
