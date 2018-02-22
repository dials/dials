from __future__ import absolute_import, division, print_function

import pytest

from cctbx import sgtbx
from scitbx import matrix
from scitbx.math import euler_angles as euler
from dxtbx.model import Crystal
from dials.algorithms.indexing import compare_orientation_matrices

def test_compare_orientation_matrices():
  # try and see if we can get back the original rotation matrix and euler angles
  real_space_a = matrix.col((10,0,0))
  real_space_b = matrix.col((0,10,10))
  real_space_c = matrix.col((0,0,10))
  euler_angles = (1.3, 5.6, 7.8)
  R = matrix.sqr(
    euler.xyz_matrix(euler_angles[0], euler_angles[1], euler_angles[2]))
  crystal_a = Crystal(real_space_a,
                      real_space_b,
                      real_space_c,
                      space_group=sgtbx.space_group('P 1'))
  crystal_b = Crystal(R * real_space_a,
                      R * real_space_b,
                      R * real_space_c,
                      space_group=sgtbx.space_group('P 1'))
  assert pytest.approx(matrix.sqr(crystal_b.get_U()) *
                       matrix.sqr(crystal_a.get_U()).transpose(), R)
  best_R_ab, best_axis, best_angle, best_cb_op = \
    compare_orientation_matrices.difference_rotation_matrix_axis_angle(
      crystal_a,
      crystal_b)
  best_euler_angles = euler.xyz_angles(best_R_ab)
  assert pytest.approx(best_euler_angles, euler_angles)
  assert best_cb_op.is_identity_op()
  assert pytest.approx(best_R_ab, R)

  # now see if we can deconvolute the original euler angles after applying
  # a change of basis to one of the crystals
  crystal_a = Crystal(real_space_a,
                      real_space_b,
                      real_space_c,
                      space_group=sgtbx.space_group('I 2 3'))
  crystal_b = Crystal(R * real_space_a,
                      R * real_space_b,
                      R * real_space_c,
                      space_group=sgtbx.space_group('I 2 3'))
  cb_op = sgtbx.change_of_basis_op('z,x,y')
  crystal_b = crystal_b.change_basis(cb_op)
  best_R_ab, best_axis, best_angle, best_cb_op = \
    compare_orientation_matrices.difference_rotation_matrix_axis_angle(
      crystal_a,
      crystal_b)
  best_euler_angles = euler.xyz_angles(best_R_ab)
  assert pytest.approx(best_euler_angles, euler_angles)
  assert best_cb_op.c() == cb_op.inverse().c()
  assert pytest.approx(best_R_ab, R)
