#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Test decomposition of rotation matrices around arbitrary axes"""

from __future__ import absolute_import, division
import math
import random
from scitbx import matrix
from libtbx.test_utils import approx_equal

from dials.algorithms.refinement.rotation_decomposition import \
  solve_r3_rotation_for_angles_given_axes

def test_rotation_matrices(phi1, phi2, phi3):

  # compose rotation matrix
  R1 = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(phi1, deg=False)
  R2 = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(phi2, deg=False)
  R3 = matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(phi3, deg=False)
  R = R1*R2*R3

  # obtain the solutions
  sol1 = solve_r3_rotation_for_angles_given_axes(R, (1,0,0), (0,1,0), (0,0,1),
    smaller_phi2_solution=False)
  sol2 = solve_r3_rotation_for_angles_given_axes(R, (1,0,0), (0,1,0), (0,0,1),
    smaller_phi2_solution=True)

  # recompose these into rotation matrices
  tstR1 = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(sol1[0],
    deg=False)
  tstR1 *= matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(sol1[1],
    deg=False)
  tstR1 *= matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(sol1[2],
    deg=False)

  tstR2 = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(sol2[0],
    deg=False)
  tstR2 *= matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(sol2[1],
    deg=False)
  tstR2 *= matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(sol2[2],
    deg=False)

  # the two solution sets must reproduce the same rotation matrix
  assert approx_equal(tstR1, tstR2)

  # additionally, this matrix must be equal to the original matrix
  assert approx_equal(tstR1, R)

  # check rotated vectors
  #vec = R * matrix.col((0,0,1))
  #tstvec1 = tstR1 * matrix.col((0,0,1))
  #tstvec2 = tstR2 * matrix.col((0,0,1))
  #print R
  #print tstR1
  #print tstR2
  #print vec
  #print tstvec1
  #print tstvec2
  return

def test_vs_euler_angles_xyz_angles(phi1, phi2, phi3):

  from scitbx.math import euler_angles_xyz_angles

  # compose rotation matrix
  R1 = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(phi1, deg=False)
  R2 = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(phi2, deg=False)
  R3 = matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(phi3, deg=False)
  R = R1*R2*R3

  # get two solution sets for the principal axes
  sol1 = solve_r3_rotation_for_angles_given_axes(R, (1,0,0), (0,1,0), (0,0,1),
    smaller_phi2_solution=False, deg=True)
  sol2 = solve_r3_rotation_for_angles_given_axes(R, (1,0,0), (0,1,0), (0,0,1),
    smaller_phi2_solution=True, deg=True)

  # get the solution set provided by euler_angles_xyz_angles
  tst = euler_angles_xyz_angles(R)

  # one of these must match
  assert any([approx_equal(sol1, tst, out=None),
              approx_equal(sol2, tst, out=None)])
  return

def random_vector():
  v = matrix.col((1,0,0)).rotate_around_origin(matrix.col((0,1,0)),
    random.uniform(0, math.pi))
  return v.rotate_around_origin(matrix.col((1,0,0)),
    random.uniform(0, 2.*math.pi)).normalize()

def test_random_axes_and_angles():

  # random axes
  e1, e2, e3 = [random_vector() for i in range(3)]

  # random angles
  phi1 = random.uniform(-math.pi, math.pi)
  phi2 = random.uniform(-math.pi, math.pi)
  phi3 = random.uniform(-math.pi, math.pi)

  # compose rotation matrix
  R1 = e1.axis_and_angle_as_r3_rotation_matrix(phi1, deg=False)
  R2 = e2.axis_and_angle_as_r3_rotation_matrix(phi2, deg=False)
  R3 = e3.axis_and_angle_as_r3_rotation_matrix(phi3, deg=False)
  R = R1*R2*R3

  # obtain solution sets
  sol1 = solve_r3_rotation_for_angles_given_axes(R, e1, e2, e3, True)
  sol2 = solve_r3_rotation_for_angles_given_axes(R, e1, e2, e3, False)

  # if either one of the solution sets is not found, indicate that this test
  # was skipped
  if sol1 is None or sol2 is None: return "skipped"

  # there are special solutions where rotations 1 and 3 are about the same axis.
  # solve_r3_rotation_for_angles_given_axes provides solutions where phi1 is
  # 0.0 and the rotation is fully expressed by phi3. Skip these degenerate
  # cases too.
  if sol1[0] == 0.0 and sol2[0] == 0.0: return "skipped"

  # one of sol1 or sol2 should match
  tst = (phi1, phi2, phi3)
  assert any([approx_equal(sol1, tst, out=None),
              approx_equal(sol2, tst, out=None)])

  return

if __name__ == '__main__':

  n_tests=100

  # tests using principal axes, and angles in the range +/-30 deg
  for i in range(n_tests):

    phi1 = random.uniform(-math.pi/6, math.pi/6)
    phi2 = random.uniform(-math.pi/6, math.pi/6)
    phi3 = random.uniform(-math.pi/6, math.pi/6)
    test_rotation_matrices(phi1, phi2, phi3)
    test_vs_euler_angles_xyz_angles(phi1, phi2, phi3)

  # tests using arbitrary axes and angles, ensuring not too many are skipped
  results = [test_random_axes_and_angles() for i in range(n_tests)]
  assert results.count("skipped") < 0.2 * n_tests

  print "OK"

