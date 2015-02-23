#!/usr/bin/env dials.python
#
# Play area trying to derive useful derivatives for Rodrigues rotation formula
# http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

sage_script = '''
r = var('r')
s = var('s')

k1 = cos(r) * sin(s)
k2 = sin(r) * sin(s)
k3 = cos(s)
K = matrix([[0, -k3, k2], [k3, 0, -k1], [-k2, k1, 0]])
t = var('t')
I = matrix.identity(3)

R = I + sin(t) * K + (1 - cos(t)) * K * K

dRdr = diff(R, r)
dRds = diff(R, s)

for i in range(3):
  for j in range(3):
    print dRdr[i,j].trig_reduce().full_simplify()

for i in range(3):
  for j in range(3):
    print dRds[i,j].trig_reduce().full_simplify()

'''

def skew_symm(v):
  '''Make matrix [v]_x from v. Essentially multiply vector by SO(3) basis
  set Lx, Ly, Lz. Equation (2) from Gallego & Yezzi paper.'''
  import scitbx.matrix

  L1 = scitbx.matrix.sqr((0, 0, 0, 0, 0, -1, 0, 1, 0))
  L2 = scitbx.matrix.sqr((0, 0, 1, 0, 0, 0, -1, 0, 0))
  L3 = scitbx.matrix.sqr((0, -1, 0, 1, 0, 0, 0, 0, 0))

  v1, v2, v3 = v.elems

  return v1 * L1 + v2 * L2 + v3 * L3

def dR_dki_gw(r, s, t):
  from scitbx import matrix
  from math import sin, cos

  cr = cos(r)
  sr = sin(r)
  cs = cos(s)
  ss = sin(s)
  ct = cos(t)
  st = sin(t)

  dRdr = matrix.sqr([
    -1/2*(cos(2*s)*sin(2*r) - sin(2*r))*ct + 1/2*cos(2*s)*sin(2*r) - cr*sr,
    cr**2 - 1/2*cos(2*r)*cos(2*s) + 1/2*(cos(2*r)*cos(2*s) - cos(2*r))*ct - 1/2,
    1/2*ct*sr*sin(2*s) + cr*ss*st - 1/2*sr*sin(2*s),
    cr**2 - 1/2*cos(2*r)*cos(2*s) + 1/2*(cos(2*r)*cos(2*s) - cos(2*r))*ct - 1/2,
    1/2*(cos(2*s)*sin(2*r) - sin(2*r))*ct - 1/2*cos(2*s)*sin(2*r) + cr*sr,
    -1/2*cr*ct*sin(2*s) + sr*ss*st + 1/2*cr*sin(2*s),
    1/2*ct*sr*sin(2*s) - cr*ss*st - 1/2*sr*sin(2*s),
    -1/2*cr*ct*sin(2*s) - sr*ss*st + 1/2*cr*sin(2*s),
    0])

  dRds = matrix.sqr([
    -1/2*(cos(2*r) + 1)*ct*sin(2*s) + 1/2*cos(2*r)*sin(2*s) + cs*ss,
    -1/2*ct*sin(2*r)*sin(2*s) + 1/2*sin(2*r)*sin(2*s) + ss*st,
    -cr*cos(2*s)*ct + cs*sr*st + cr*cos(2*s),
    -1/2*ct*sin(2*r)*sin(2*s) + 1/2*sin(2*r)*sin(2*s) - ss*st,
    1/2*(cos(2*r) - 1)*ct*sin(2*s) - 1/2*cos(2*r)*sin(2*s) + cs*ss,
    -cos(2*s)*ct*sr - cr*cs*st + cos(2*s)*sr,
    -cr*cos(2*s)*ct - cs*sr*st + cr*cos(2*s),
    -cos(2*s)*ct*sr + cr*cs*st + cos(2*s)*sr,
    ct*sin(2*s) - 2*cs*ss])

  return dRdr, dRds

def Rtk(t, k):
  import scitbx.matrix
  import math

  I = scitbx.matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  c = math.cos(t)
  s = math.sin(t)
  K = skew_symm(k)

  R = I + s * K + (1 - c) * K * K

  # make sure answer is right
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

def dx():
  import random
  import math
  import scitbx.matrix

  # pick random rotation axis & angle - pick spherical coordinate basis to make
  # life easier in performing finite difference calculations

  t = math.pi * random.random()
  r = 2 * math.pi * random.random()
  s = math.pi * random.random()
  k = krs(r, s)

  # perform finite difference calculation of dR/dk1 as
  # (1/2 eps) * (R(k1+eps) - R(k1-eps))

  eps = 1.0e-8

  Rpls = Rtk(t, krs(r + eps, s))
  Rmns = Rtk(t, krs(r - eps, s))
  ndRr = (0.5 / eps) * (Rpls - Rmns)

  Rpls = Rtk(t, krs(r, s + eps))
  Rmns = Rtk(t, krs(r, s - eps))
  ndRs = (0.5 / eps) * (Rpls - Rmns)

  # now call on derivative calculation
  dRr, dRs = dR_dki_gw(r, s, t)

  print 'Numerical dR/dr:'
  print ndRr

  print 'Sums dR/dr:'
  print dRr

  print 'Numerical dR/ds:'
  print ndRs

  print 'Sums dR/ds:'
  print dRs

  I = scitbx.matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  print 'Summed differences:'
  print sum((ndRr - dRr).elems)
  print sum((ndRs - dRs).elems)

dx()
