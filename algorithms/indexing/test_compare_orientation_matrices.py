from __future__ import absolute_import, division, print_function

import pytest
from cctbx import sgtbx
from scitbx import matrix
from scitbx.math import euler_angles as euler
from dxtbx.model import Crystal
from dials.algorithms.indexing import compare_orientation_matrices


def test_compare_orientation_matrices():
    # try and see if we can get back the original rotation matrix and euler angles
    real_space_a = matrix.col((10, 0, 0))
    real_space_b = matrix.col((0, 10, 10))
    real_space_c = matrix.col((0, 0, 10))
    euler_angles = (1.3, 5.6, 7.8)
    R = matrix.sqr(euler.xyz_matrix(euler_angles[0], euler_angles[1], euler_angles[2]))
    crystal_a = Crystal(
        real_space_a, real_space_b, real_space_c, space_group=sgtbx.space_group("P 1")
    )
    crystal_b = Crystal(
        R * real_space_a,
        R * real_space_b,
        R * real_space_c,
        space_group=sgtbx.space_group("P 1"),
    )
    assert (
        matrix.sqr(crystal_b.get_U()) * matrix.sqr(crystal_a.get_U()).transpose()
    ).elems == pytest.approx(R.elems)
    best_R_ab, best_axis, best_angle, best_cb_op = compare_orientation_matrices.difference_rotation_matrix_axis_angle(
        crystal_a, crystal_b
    )
    best_euler_angles = euler.xyz_angles(best_R_ab)
    assert best_euler_angles == pytest.approx(euler_angles)
    assert best_cb_op.is_identity_op()
    assert best_R_ab.elems == pytest.approx(R.elems)

    # now see if we can deconvolute the original euler angles after applying
    # a change of basis to one of the crystals
    crystal_a = Crystal(
        real_space_a, real_space_b, real_space_c, space_group=sgtbx.space_group("I 2 3")
    )
    cb_op = sgtbx.change_of_basis_op("z,x,y")
    crystal_b = Crystal(
        R * real_space_a,
        R * real_space_b,
        R * real_space_c,
        space_group=sgtbx.space_group("I 2 3"),
    ).change_basis(cb_op)
    best_R_ab, best_axis, best_angle, best_cb_op = compare_orientation_matrices.difference_rotation_matrix_axis_angle(
        crystal_a, crystal_b
    )
    best_euler_angles = euler.xyz_angles(best_R_ab)
    assert best_euler_angles == pytest.approx(euler_angles)
    assert best_cb_op.c() == cb_op.inverse().c()
    assert best_R_ab.elems == pytest.approx(R.elems)

    crystal_c = crystal_b.change_basis(sgtbx.change_of_basis_op("-y,-z,x"))
    assert crystal_c != crystal_b

    s = compare_orientation_matrices.rotation_matrix_differences(
        [crystal_a, crystal_b, crystal_c], comparison="pairwise"
    )
    s = "\n".join(s.splitlines()[:-1]).replace("-0.000", "0.000")
    print(s)
    assert (
        s
        == """\
Change of basis op: b,c,a
Rotation matrix to transform crystal 1 to crystal 2:
{{0.986, -0.135, 0.098},
 {0.138, 0.990, -0.023},
 {-0.094, 0.036, 0.995}}
Rotation of 9.738 degrees about axis (0.172, 0.565, 0.807)

Change of basis op: -a,-b,c
Rotation matrix to transform crystal 1 to crystal 3:
{{0.986, -0.135, 0.098},
 {0.138, 0.990, -0.023},
 {-0.094, 0.036, 0.995}}
Rotation of 9.738 degrees about axis (0.172, 0.565, 0.807)

Change of basis op: c,-a,-b
Rotation matrix to transform crystal 2 to crystal 3:
{{1.000, 0.000, 0.000},
 {0.000, 1.000, 0.000},
 {0.000, 0.000, 1.000}}"""
    )

    s = compare_orientation_matrices.rotation_matrix_differences(
        [crystal_a, crystal_b, crystal_c], comparison="sequential"
    )
    s = "\n".join(s.splitlines()[:-1]).replace("-0.000", "0.000")
    print(s)
    assert (
        s
        == """\
Change of basis op: b,c,a
Rotation matrix to transform crystal 1 to crystal 2:
{{0.986, -0.135, 0.098},
 {0.138, 0.990, -0.023},
 {-0.094, 0.036, 0.995}}
Rotation of 9.738 degrees about axis (0.172, 0.565, 0.807)

Change of basis op: c,-a,-b
Rotation matrix to transform crystal 2 to crystal 3:
{{1.000, 0.000, 0.000},
 {0.000, 1.000, 0.000},
 {0.000, 0.000, 1.000}}"""
    )

    s = compare_orientation_matrices.rotation_matrix_differences(
        (crystal_a, crystal_b), miller_indices=((1, 0, 0), (1, 1, 0))
    )
    assert (
        s
        == """\
Change of basis op: b,c,a
Rotation matrix to transform crystal 1 to crystal 2:
{{0.986, -0.135, 0.098},
 {0.138, 0.990, -0.023},
 {-0.094, 0.036, 0.995}}
Rotation of 9.738 degrees about axis (0.172, 0.565, 0.807)
(1,0,0): 15.26 deg
(1,1,0): 9.12 deg
"""
    )
