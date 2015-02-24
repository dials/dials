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

dx()
