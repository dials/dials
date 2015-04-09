#!/usr/bin/env python
#
#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import math
from scitbx import matrix

"""A Python version of the Rotation::euler_explicit method from Pointless"""

def solve_r3_rotation_for_angles_given_axes(R, e1, e2, e3,
  smaller_phi2_solution=True, deg=False):
  """decompose 3*3 rotation matrix R into three rotations about the axes
  e1, e2 and e3

  return angles given axes; return None if no solution
  smaller_phi2_solution = True or False for the two solutions

  Similar to scitbx.math.euler_angles, except the axes here are arbitrary.

  If the angle between e1 & e2 or e2 & e3 are < 90 degrees, then
  there are matrices which have no solution

  The algorithm for decomposing a rotation into three arbitrary rotations was
  described by David Thomas (Thomas, D.J.(1990) Acta Cryst. A46, 321-343)
  and detailed by Gerard Bricogne (Proceedings of the CCP4 Study Weekend 23-24
  January 1987).

  This implementation is a translation of the C++ code Rotation::euler_explicit
  from Pointless and Aimless by Phil Evans.

  The notation used follows Bricogne

  R =  R(e1, phi1) R(e2, phi2) R(e3, phi3)

  where R(e, phi) is the rotation matrix for a rotation by angle phi around
  axis e"""

  assert R.is_r3_rotation_matrix()
  e1 = matrix.col(e1).normalize()
  e2 = matrix.col(e2).normalize()
  e3 = matrix.col(e3).normalize()
  # Fail if e2 & e3 are parallel
  e2xe3 = e2.cross(e3)
  if e2xe3.length_sq() < 1.0e-6: return None
  # Make a unit test vector
  u = e2xe3.normalize()
  e1e2 = e1.dot(e2)
  e1e3 = e1.dot(e3)
  e2e3 = e2.dot(e3)
  e1e2e3 = e1.dot(e2.cross(e3))
  Re3 = R*e3
  e1Re3 = e1.dot(Re3)
  # ** Step 1 ** Calculation of phi2  (Bricogne equation (4))
  # e1.(R e3) = (e1.e2)(e2.e3) + {(e1.e3) - (e1.e2)(e2.e3)} cos(phi2)
  #           + (e1.e2 x e3) sin(phi2)
  # The coefficients of cos & sin phi2
  cc = e1e3 - e1e2 * e2e3
  ss = e1e2e3
  # Fail if both are zero (indeterminate)
  if abs(cc) < 1.0e-6 and abs(ss) < 1.0e-6: return None
  norm = math.sqrt(cc*cc + ss*ss)
  rhs = (e1Re3 - e1e2 * e2e3) / norm
  # abs(rhs) should not be greater than 1.0, allowing a small tolerance
  if abs(rhs) > 1.000002: return None
  if rhs > 1.0: rhs = 1.0
  elif rhs < -1.0: rhs = -1.0
  cc /= norm
  ss /= norm
  # Solve rhs = cos(phi2) * cc + sin(phi2) * ss
  # using cos(a-b) = cos(a) cos(b) + sin(a) sin(b)
  # where b = phi2
  a = math.atan2(ss, cc)
  amb = math.acos(rhs)
  # Two solutions in range -pi to +pi
  # Note that if e1 == e3, ss = 0, a = 0 & phi2b = -phi2a
  phi2a = a - amb
  if phi2a > math.pi: phi2a -= 2.0*math.pi
  elif phi2a < -math.pi: phi2a += 2.0*math.pi
  phi2b = a + amb
  if phi2b > math.pi: phi2b -= 2.0*math.pi
  elif phi2b < -math.pi: phi2b += 2.0*math.pi
  if smaller_phi2_solution:
    if abs(phi2a) < abs(phi2b):
      phi2 = phi2a
    else:
      phi2 = phi2b
  else:
    if abs(phi2a) > abs(phi2b):
      phi2 = phi2a
    else:
      phi2 = phi2b
  # ** Step 2 ** Calculation of phi1
  R2 = e2.axis_and_angle_as_r3_rotation_matrix(phi2, deg=False)
  R2inv = R2.transpose()
  v = R2 * e3
  w = Re3
  v1 = v - (v.dot(e1)) * e1
  w1 = w - (w.dot(e1)) * e1
  norm = v1.dot(v1)*w1.dot(w1)
  # If norm = 0, rotations 1 & 3 are around same axis (for this phi2),
  # so any value for phi1 is OK
  if (norm > 1.0e-8):
    norm = math.sqrt(norm)
    # cos(phi1) = (v1.w1)/norm
    # sin(phi1) = (v1.w1 x e1)/norm
    phi1 = math.atan2(v1.dot(w1.cross(e1))/norm, v1.dot(w1)/norm)
    if phi1 >  math.pi: phi1 -= 2.0*math.pi
    if phi1 < -math.pi: phi1 += 2.0*math.pi
  else:
    phi1 = 0.0
  # ** Step 3 ** Calculation of phi3
  R1inv = e1.axis_and_angle_as_r3_rotation_matrix(-1.*phi1, deg=False)
  R3 = R2inv * R1inv * R
  R3u = R3 * u;
  # sin(phi3) = u.R3u x e3
  # cos(phi3) = u.R3u
  phi3 = math.atan2(u.dot(R3u.cross(e3)), u.dot(R3u))
  if deg:
    phi1, phi2, phi3 = tuple([x * 180/math.pi for x in (phi1, phi2, phi3)])
  return phi1, phi2, phi3
