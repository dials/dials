from __future__ import division
from libtbx.test_utils import approx_equal

def exercise():
  from dials.algorithms.indexing import compare_orientation_matrices
  from dxtbx.model.crystal import crystal_model
  from cctbx import sgtbx
  from scitbx import matrix
  from scitbx.math import euler_angles as euler

  # try and see if we can get back the original rotation matrix and euler angles
  real_space_a = matrix.col((10,0,0))
  real_space_b = matrix.col((0,10,10))
  real_space_c = matrix.col((0,0,10))
  euler_angles = (1.3, 5.6, 7.8)
  R = matrix.sqr(
    euler.xyz_matrix(euler_angles[0], euler_angles[1], euler_angles[2]))
  crystal_a = crystal_model(real_space_a,
                            real_space_b,
                            real_space_c,
                            space_group=sgtbx.space_group('P 1'))
  crystal_b = crystal_model(R * real_space_a,
                            R * real_space_b,
                            R * real_space_c,
                            space_group=sgtbx.space_group('P 1'))
  assert approx_equal(crystal_b.get_U() * crystal_a.get_U().transpose(), R)
  best_R_ab, best_axis, best_angle, best_cb_op = \
    compare_orientation_matrices.difference_rotation_matrix_axis_angle(
      crystal_a,
      crystal_b)
  best_euler_angles = euler.xyz_angles(best_R_ab)
  assert approx_equal(best_euler_angles, euler_angles)
  assert best_cb_op.is_identity_op()
  assert approx_equal(best_R_ab, R)

  # now see if we can deconvolute the original euler angles after applying
  # a change of basis to one of the crystals
  crystal_a = crystal_model(real_space_a,
                            real_space_b,
                            real_space_c,
                            space_group=sgtbx.space_group('I 2 3'))
  crystal_b = crystal_model(R * real_space_a,
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
  assert approx_equal(best_euler_angles, euler_angles)
  assert best_cb_op.c() == cb_op.inverse().c()
  assert approx_equal(best_R_ab, R)




def run():
  exercise()
  print "OK"

if __name__ == '__main__':
  run()
