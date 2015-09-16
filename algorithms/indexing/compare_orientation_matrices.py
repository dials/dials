from __future__ import division
from scitbx.math import euler_angles as euler
from scitbx import matrix


def difference_rotation_matrix_and_euler_angles(crystal_a, crystal_b):
  from cctbx import sgtbx
  #assert crystal_a.get_space_group() == crystal_b.get_space_group()
  space_group = crystal_b.get_space_group()
  min_sum_euler_angles = 1e8
  best_R_ab = None
  best_cb_op = None
  # iterate over space group ops to find smallest differences
  for i_op, op in enumerate(space_group.build_derived_laue_group().all_ops()):
    if op.r().determinant() < 0:
      continue
    elif not op.t().is_zero():
      continue
    cb_op = sgtbx.change_of_basis_op(op.inverse())
    crystal_b_sym = crystal_b.change_basis(cb_op)
    U_a = crystal_a.get_U()
    U_b = crystal_b_sym.get_U()
    assert U_a.is_r3_rotation_matrix()
    assert U_b.is_r3_rotation_matrix()
    # the rotation matrix to transform from U_a to U_b
    R_ab = U_b * U_a.transpose()
    euler_angles = euler.xyz_angles(R_ab)
    sum_euler_angles = sum(abs(ea) for ea in euler_angles)
    if sum_euler_angles < min_sum_euler_angles:
      min_sum_euler_angles = sum_euler_angles
      best_R_ab = R_ab
      best_cb_op = cb_op

  best_euler_angles = euler.xyz_angles(best_R_ab)
  return best_R_ab, best_euler_angles, best_cb_op


def show_rotation_matrix_differences(crystal_models, out=None, miller_indices=None):
  if out is None:
    import sys
    out = sys.stdout
  for i in range(len(crystal_models)):
    for j in range(i+1, len(crystal_models)):
      R_ij, euler_angles, cb_op = difference_rotation_matrix_and_euler_angles(
        crystal_models[i], crystal_models[j])
      print "Rotation matrix to transform crystal %i to crystal %i" %(
        i+1, j+1)
      print "Change of basis op: %s" %cb_op
      print R_ij.mathematica_form(format="%.3f", one_row_per_line=True)
      print "Euler angles (xyz): %.2f, %.2f, %.2f" %euler_angles
      if miller_indices is not None:
        for hkl in miller_indices:
          cm_i = crystal_models[i]
          cm_j = crystal_models[j].change_basis(cb_op)
          A_i = cm_i.get_A()
          A_j = cm_j.get_A()
          a_star_i = matrix.col(A_i[0:3])
          b_star_i = matrix.col(A_i[3:6])
          c_star_i = matrix.col(A_i[6:9])
          a_star_j = matrix.col(A_j[0:3])
          b_star_j = matrix.col(A_j[3:6])
          c_star_j = matrix.col(A_j[6:9])

          v_i = hkl[0] * a_star_i + hkl[1] * b_star_i + hkl[2] * c_star_i
          v_j = hkl[0] * a_star_j + hkl[1] * b_star_j + hkl[2] * c_star_j

          print '(%i,%i,%i):' %hkl, '%.2f deg' %v_i.angle(v_j, deg=True)
      print
