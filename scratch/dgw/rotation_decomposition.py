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
from math import sqrt, acos, atan2, pi
from scitbx import matrix

"""Python conversion of Rotation::euler_explicit method from Pointless"""
TWOPI = 2.0*pi

def angle_in_range(angle):
  # Put angle in range -pi to +pi
  if angle >  pi: return angle - TWOPI
  if angle < -pi: return angle + TWOPI
  return angle;

def decompose_to_angles_given_axes(R, e1, e2, e3, solution=1):
  """decompose 3*3 rotation matrix R into three rotations about the axes
  e1, e2 and e3

  return angles given axes; return None if no solution
  solution = 1 or 2 for the two solutions

  If the angle between e1 & e2 or e2 & e3 are < 90 degrees, then
  there are matrices which have no solution

  Note that the algorithm for decomposing a rotation into three arbitrary
  rotations was described by David Thomas
  (Thomas, D.J.(1990) Acta Cryst. A46, 321-343)
  and detailed by Gerard Bricogne
  (Proceedings of the CCP4 Study Weekend 23-24 January 1987)
  The code here is a essentially translation of a Fortran routine
  by Zbyszek Otwinowski  (solveu)

  The notation used follows Bricogne

   R =  R(e1, phi1) R(e2, phi2) R(e3, phi3)

   where R(e, phi) is the rotation matrix for a rotation by angle
  phi around axis e"""

  assert R.is_r3_rotation_matrix()
  e1 = matrix.col(e1)
  e2 = matrix.col(e2)
  e3 = matrix.col(e3)

  e1e2 = e1.dot(e2)
  e1e3 = e1.dot(e3)
  e2e3 = e1.dot(e3)
  e1e2e3 = e1.dot(e2.cross(e3))

  # R.e3
  Re3 = R*e3
  # e1. Re3
  e1Re3 = e1.dot(Re3)

  # We need a test vector perpendicular to e3: use e2 x e3
  u = e2.cross(e3)
  # Fail if e2 & e3 are parallel
  if u.length_sq() < 1.0e-6: return None

  # ** Step 1 ** Calculation of phi2  (Bricogne equation (4))
  # e1.(R e3) = (e1.e2)(e2.e3) + {(e1.e3) - (e1.e2)(e2.e3)} cos(phi2)
  #           + (e1.e2 x e3) sin(phi2)
  # The coefficients of cos & sin phi2
  cc = e1e3 - e1e2 * e2e3
  ss = e1e2e3
  # Fail if both are zero (indeterminate)
  if abs(cc) < 1.0e-6 and abs(ss) < 1.0e-6: return None
  # Normalise equation (4)
  norm = sqrt(cc*cc + ss*ss)
  # rhs = e1.Re3 - (e1.e2)(e2.e3)
  rhs = (e1Re3 - e1e2 * e2e3) / norm
  # abs(rhs) should not be greater than 1.0
  if abs(rhs) > 1.000002: return None
  # Allow a small tolerance
  if rhs > 1.0: rhs = 1.0
  elif rhs < -1.0: rhs = -1.0
  cc /= norm
  ss /= norm
  # Solve  rhs  = cos(phi2) * cc + sin(phi2) * ss
  # using cos(a-b) = cos(a) cos(b) + sin(a) sin(b)
  # where b = phi2
  a = atan2(ss, cc)
  amb = acos(rhs)  # +-(a-b)
  # Two solutions (from arc cos), phi2 = b = a -+(a-b)
  #   1)  phi2 = a - amb
  #   2)  phi2 = a + amb
  # in range -pi to +pi
  # Note that if e1 == e3, ss = 0, a = 0 & phi2b = -phi2a
  phi2a = angle_in_range(a - amb)
  phi2b = angle_in_range(a + amb)
  # Choose solution 1 (larger phi2 ie positive) or
  # solution 2 (smaller (negative) phi2)
  #  make phi2a > phi2b
  if phi2a < phi2b: phi2a, phi2b = phi2b, phi2a
  phi2 = phi2b
  if (solution == 1): phi2 = phi2a;   # pick one solution

  # ** Step 2 ** Calculation of phi1
  phi1 = 0.0;
  R2 = e2.axis_and_angle_as_r3_rotation_matrix(phi2, deg=False)
  R2inv = e2.axis_and_angle_as_r3_rotation_matrix(-1.0*phi2, deg=False)
  v = R2 * e3
  w = Re3
  # v1 = v - (v.e1) e1
  # w1 = w - (w.e1) e1
  v1 = v - (v.dot(e1)) * e1
  w1 = w - (w.dot(e1)) * e1
  norm = v1.dot(v1)*w1.dot(w1)
  # If norm = 0, rotations 1 & 3 are around same axis (for this phi2),
  # so any value for phi1 is OK (leave = 0.0)
  if (norm > 1.0e-8):
    norm = sqrt(norm)
    # norm = sqrt((v1.v1)(w1.w1))
    # cos(phi1) = (v1.w1)/norm
    # sin(phi1) = (v1.w1 x e1)/norm
    #phi1 = angle_in_range(atan2((v1*Vec3<ftype>::cross(w1, e1))/norm, (v1*w1)/norm));
    phi1 = angle_in_range(atan2(v1.dot(w1.cross(e1))/norm, v1.dot(w1)/norm))
  R1inv = e1.axis_and_angle_as_r3_rotation_matrix(-1.*phi1, deg=False)

  # ** Step 3 ** Calculation of phi3
  R3 = R2inv * R1inv * R
  R3u = R3 * u;
  #  sin(phi3) = u. R3u x e3
  #  cos(phi3) = u . R3u
  phi3 = atan2(u.dot(R3u.cross(e3)), u.dot(R3u))
  return phi1, phi2, phi3;

if __name__ == '__main__':

  R1=matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(pi/7, deg=False)
  R2=matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(pi/9, deg=False)
  R3=matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(-pi/8, deg=False)

  print "target angles"
  print pi/7, pi/9, -pi/8
  R=R1*R2*R3

  print "solution 1"
  print decompose_to_angles_given_axes(R, (1,0,0), (0,1,0), (0,0,1), 1)
  print "solution 2"
  print decompose_to_angles_given_axes(R, (1,0,0), (0,1,0), (0,0,1), 2)
