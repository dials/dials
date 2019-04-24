from __future__ import absolute_import, division, print_function

import math

from scitbx import matrix
from scitbx.math import r3_rotation_axis_and_angle_from_matrix


def difference_rotation_matrix_axis_angle(crystal_a, crystal_b, target_angle=0):
    from cctbx import sgtbx

    # assert crystal_a.get_space_group() == crystal_b.get_space_group()
    space_group = crystal_b.get_space_group()
    best_R_ab = None
    best_cb_op = None
    best_axis = None
    best_angle = 1e8
    # iterate over space group ops to find smallest differences
    for i_op, op in enumerate(space_group.build_derived_laue_group().all_ops()):
        if op.r().determinant() < 0:
            continue
        elif not op.t().is_zero():
            continue
        cb_op = sgtbx.change_of_basis_op(op.inverse())
        crystal_b_sym = crystal_b.change_basis(cb_op)
        U_a = matrix.sqr(crystal_a.get_U())
        U_b = matrix.sqr(crystal_b_sym.get_U())
        assert U_a.is_r3_rotation_matrix()
        assert U_b.is_r3_rotation_matrix()
        # the rotation matrix to transform from U_a to U_b
        R_ab = U_b * U_a.transpose()
        axis_angle = r3_rotation_axis_and_angle_from_matrix(R_ab)
        axis = axis_angle.axis
        angle = axis_angle.angle() * 180.0 / math.pi
        for sign in (+1, -1):
            if abs(sign * angle - target_angle) < abs(best_angle - target_angle):
                best_angle = sign * angle
                best_axis = tuple(sign * a for a in axis)
                best_R_ab = R_ab if sign > 0 else R_ab.inverse()
                best_cb_op = cb_op if sign > 0 else cb_op.inverse()

    return best_R_ab, best_axis, best_angle, best_cb_op


def rotation_matrix_differences(
    crystal_models, miller_indices=None, comparison="pairwise"
):
    assert comparison in ("pairwise", "sequential"), comparison
    output = []
    for i in range(len(crystal_models)):
        for j in range(i + 1, len(crystal_models)):
            if comparison == "sequential" and j > i + 1:
                break
            R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
                crystal_models[i], crystal_models[j]
            )
            output.append("Change of basis op: %s" % cb_op)
            output.append(
                "Rotation matrix to transform crystal %i to crystal %i:"
                % (i + 1, j + 1)
            )
            output.append(R_ij.mathematica_form(format="%.3f", one_row_per_line=True))
            output.append(
                ("Rotation of %.3f degrees " % angle)
                + ("about axis (%.3f, %.3f, %.3f)" % axis)
            )
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

                    output.append(
                        ("(%i,%i,%i): " % hkl) + ("%.2f deg" % v_i.angle(v_j, deg=True))
                    )
            output.append("")
    return "\n".join(output)
