def skew_symm(v):
  '''Make matrix [v]_x from v. Essentially multiply vector by SO(3) basis
  set Lx, Ly, Lz. Equation (2) from Gallego & Yezzi paper.'''
  import scitbx.matrix

  L1 = scitbx.matrix.sqr((0, 0, 0, 0, 0, -1, 0, 1, 0))
  L2 = scitbx.matrix.sqr((0, 0, 1, 0, 0, 0, -1, 0, 0))
  L3 = scitbx.matrix.sqr((0, -1, 0, 1, 0, 0, 0, 0, 0))

  v1, v2, v3 = v.elems

  return v1 * L1 + v2 * L2 + v3 * L3

def dR_dki_7(t, k):
  '''Compute derivative matrices dR/dk_i for rotation about axis k of
  angle t. Equation (7) from Gallego & Yezzi paper, though that is derivative
  of R as function of v = t * k => remove one factor of t = |v|.'''

  import scitbx.matrix
  import math

  K = skew_symm(k)
  I = scitbx.matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  e1 = scitbx.matrix.col((1, 0, 0))
  e2 = scitbx.matrix.col((0, 1, 0))
  e3 = scitbx.matrix.col((0, 0, 1))

  e = (e1, e2, e3)

  c = math.cos(t)
  s = math.sin(t)

  result = []

  for j in range(3):
    result.append((t * c * k.elems[j] * K + t * s * k.elems[j] * K * K +
                   s * skew_symm(e[j] - k.elems[j] * k) +
                   (1 - c) * (e[j] * k.transpose() +
                              k * e[j].transpose() -
                              2 * k.elems[j] * k * k.transpose())))
  return result

def dR_dki_9(t, k):
  '''Compute derivative matrices dR/dk_i for rotation about axis k of
  angle t. Equation (9) from Gallego & Yezzi paper, though that is derivative
  of R as function of v = t * k => remove one factor of t = |v|.'''

  import scitbx.matrix

  R = Rtk(t, k)
  v = t * k
  V = skew_symm(v)
  I = scitbx.matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  e1 = scitbx.matrix.col((1, 0, 0))
  e2 = scitbx.matrix.col((0, 1, 0))
  e3 = scitbx.matrix.col((0, 0, 1))

  e = (e1, e2, e3)

  result = []

  for j in range(3):
    result.append((v.elems[j] * V + skew_symm(v.cross((I - R) * e[j]))) * R / t)

  return result

def Rtk(t, k):
  import scitbx.matrix
  import math

  I = scitbx.matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  c = math.cos(t)
  s = math.sin(t)
  K = skew_symm(k)

  R = I + s * K + (1 - c) * K * K

  # make sure answer is right

  zero = R.inverse() * k.axis_and_angle_as_r3_rotation_matrix(t) - I
  assert sum(zero.elems) < 1.0e-8

  return R

def dx():
  import random
  import math
  import scitbx.matrix

  # pick random rotation axis & angle - does around the houses to get values
  # to use in finite difference calculation

  t = math.pi * random.random()
  k = scitbx.matrix.col((random.random(),
                         random.random(),
                         random.random())).normalize()
  k1, k2, k3 = k.elems

  # perform finite difference calculation of dR/dk1 as
  # (1/2 eps) * (R(k1+eps) - R(k1-eps))

  eps = 1.0e-8

  Rpls = Rtk(t, scitbx.matrix.col((k1 + eps, k2, k3)).normalize())
  Rmns = Rtk(t, scitbx.matrix.col((k1 - eps, k2, k3)).normalize())

  ndR1 = (0.5 / eps) * (Rpls - Rmns)

  print 'Numerical:'
  print ndR1

  # now call on derivative calculation
  dR1_7 = dR_dki_7(t, k)[0]
  dR1_9 = dR_dki_9(t, k)[0]

  print dR1_7 * dR1_9.inverse()


  print 'Sums:'
  print dR1_9

dx()
