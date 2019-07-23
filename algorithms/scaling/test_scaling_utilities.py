"""
Tests for scaling utilities module.
"""
from __future__ import absolute_import, division, print_function

from math import sqrt, pi
import os

import numpy as np
import pytest
from dxtbx.model import Experiment, Crystal
from dials.array_family import flex
from dials.algorithms.scaling.scaling_utilities import (
    calc_crystal_frame_vectors,
    calc_theta_phi,
    align_rotation_axis_along_z,
    set_wilson_outliers,
    quasi_normalisation,
    calculate_prescaling_correction,
    Reasons,
)
from dials_scaling_ext import (
    calculate_harmonic_tables_from_selections,
    create_sph_harm_lookup_table,
    create_sph_harm_table,
    calc_lookup_index,
)
from libtbx.test_utils import approx_equal
from mock import Mock


@pytest.fixture(scope="module")
def mock_exp():
    """Create a mock experiments object."""
    exp = Mock()
    exp.beam.get_s0.return_value = (1.0, 0.0, 0.0)
    exp.goniometer.get_rotation_axis.return_value = (0.0, 0.0, 1.0)
    return exp


@pytest.fixture(scope="module")
def test_exp_E2():
    """Create a mock experiments object."""
    exp = Experiment()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 2.0],
        "space_group_hall_symbol": " P 1",
    }
    crystal = Crystal.from_dict(exp_dict)
    exp.crystal = crystal
    return exp


@pytest.fixture(scope="module")
def test_exp_P1():
    """Create a mock experiments object."""
    exp = Experiment()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 1.0],
        "space_group_hall_symbol": " P 1",
    }
    crystal = Crystal.from_dict(exp_dict)
    exp.crystal = crystal
    return exp


@pytest.fixture(scope="module")
def test_reflection_table():
    """Return a test reflection table."""
    return generate_reflection_table()


@pytest.fixture(scope="module")
def wilson_test_reflection_table():
    """Return a test reflection table."""
    rt = flex.reflection_table()
    rt["centric_flag"] = flex.bool([True, True, False, False])
    rt["Esq"] = flex.double([50.0, 10.0, 50.0, 10.0])
    return rt


def generate_reflection_table():
    """Create a reflection table with s1 and phi."""
    rt = flex.reflection_table()
    s1_vec = (1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0))
    rt["s1"] = flex.vec3_double([s1_vec, s1_vec, s1_vec])
    rt["phi"] = flex.double([0.0, 45.0, 90.0])
    return rt


@pytest.fixture
def simple_reflection_table():
    """Create a small reflection table"""
    refl = flex.reflection_table()
    refl["intensity"] = flex.double([1.0, 2.0, 3.0])
    refl["d"] = flex.double([1.0, 2.0, 3.0])
    refl["miller_index"] = flex.miller_index([(0, 0, 3), (0, 0, 2), (0, 0, 1)])
    refl.set_flags(flex.bool(refl.size(), False), refl.flags.bad_for_scaling)
    return refl


def refl_for_norm():
    """Create 11000 refelctions in 10 groups of 1100 approx equally spaced in
    resolution."""
    intensity_array = flex.double([])
    miller_indices = flex.miller_index([])
    # a set of miller indices with h2 + k2 + l2 = [2,3,4,5,6,8,9,10,11,12],
    # which should split nicely into 10 resolution groups.
    miller_array_list = [
        (1, 1, 0),
        (1, 1, 1),
        (2, 0, 0),
        (2, 1, 0),
        (2, 1, 1),
        (2, 2, 0),
        (2, 2, 1),
        (3, 0, 1),
        (3, 1, 1),
        (2, 2, 2),
    ]
    for i in range(1, 11):
        miller_indices.extend(flex.miller_index(1100, miller_array_list[i - 1]))
        intensity_array.extend(
            flex.double(np.linspace(90, 110, num=1100, endpoint=True))
        )
    rt = flex.reflection_table()
    rt["intensity"] = intensity_array
    rt["miller_index"] = miller_indices
    rt.set_flags(flex.bool(11000, False), rt.flags.bad_for_scaling)
    return rt


def test_quasi_normalisation(simple_reflection_table, test_exp_E2, test_exp_P1):
    """Test the quasi_normalisation function."""
    # Test that for small datasets, all Esq values are set to one.
    refl = quasi_normalisation(simple_reflection_table, test_exp_E2)
    assert list(refl["Esq"]) == [1.0, 1.0, 1.0]

    rt = refl_for_norm()
    new_rt = quasi_normalisation(rt, test_exp_P1)
    for i in range(0, 9):
        assert list(new_rt["Esq"][i * 1100 : (i + 1) * 1100]) == pytest.approx(
            list(np.linspace(0.9, 1.1, num=1100, endpoint=True))
        )
    # FIXME Note, this test should actually be for i in range(0, 10):, however
    # the binner appears to create an extra bin above the highest data,
    # and then the call to interpolate causes the last values to be incorrect.


def test_calc_crystal_frame_vectors(test_reflection_table, mock_exp):
    """Test the namesake function, to check that the vectors are correctly rotated
    into the crystal frame."""
    rt, exp = test_reflection_table, mock_exp
    s0_vec = (1.0, 0.0, 0.0)
    s1_vec = (1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0))
    reflection_table = calc_crystal_frame_vectors(rt, exp)
    assert list(reflection_table["s0"]) == list(
        flex.vec3_double([s0_vec, s0_vec, s0_vec])
    )
    assert approx_equal(
        list(reflection_table["s0c"]),
        list(
            flex.vec3_double(
                [s0_vec, (1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0.0), (0.0, -1.0, 0.0)]
            )
        ),
    )
    assert approx_equal(
        list(reflection_table["s1c"]),
        list(
            flex.vec3_double(
                [
                    s1_vec,
                    (1.0 / 2.0, -1.0 / 2.0, 1.0 / sqrt(2.0)),
                    (0.0, -1.0 / sqrt(2.0), 1.0 / sqrt(2.0)),
                ]
            )
        ),
    )


def test_align_rotation_axis_along_z():
    """Test the function to rotate the coordinate system such that the rotation
    axis is along z. In the test, the rotation axis is x, so we expect the
    transformation to be: x > z, y > y, z > -x, x+z > -x+z."""
    rot_axis = flex.vec3_double([(1.0, 0.0, 0.0)])
    vectors = flex.vec3_double(
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)]
    )
    rotated_vectors = align_rotation_axis_along_z(rot_axis, vectors)
    assert approx_equal(
        list(rotated_vectors),
        list(
            flex.vec3_double(
                [(0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (-1.0, 0.0, 1.0)]
            )
        ),
    )


def test_create_sph_harm_table(test_reflection_table, mock_exp):
    """Simple test for the spherical harmonic table, constructing the table step
    by step, and verifying the values of a few easy-to-calculate entries.
    This also acts as a test for the calc_theta_phi function as well."""

    rt, exp = test_reflection_table, mock_exp
    reflection_table = calc_crystal_frame_vectors(rt, exp)
    theta_phi = calc_theta_phi(reflection_table["s0c"])
    assert approx_equal(
        list(theta_phi),
        [(pi / 2.0, 0.0), (pi / 2.0, -1.0 * pi / 4.0), (pi / 2.0, -1.0 * pi / 2.0)],
    )
    theta_phi_2 = calc_theta_phi(reflection_table["s1c"])
    assert approx_equal(
        list(theta_phi_2),
        [(pi / 4.0, 0.0), (pi / 4.0, -1.0 * pi / 4.0), (pi / 4.0, -1.0 * pi / 2.0)],
    )
    sph_h_t = create_sph_harm_table(theta_phi, theta_phi_2, 2)
    Y10 = ((3.0 / (8.0 * pi)) ** 0.5) / 2.0
    Y20 = -1.0 * ((5.0 / (256.0 * pi)) ** 0.5)
    assert approx_equal(sph_h_t[1, 0], Y10)
    assert approx_equal(sph_h_t[1, 1], Y10)
    assert approx_equal(sph_h_t[1, 2], Y10)
    assert approx_equal(sph_h_t[5, 0], Y20)
    assert approx_equal(sph_h_t[5, 1], Y20)
    assert approx_equal(sph_h_t[5, 2], Y20)
    # Now test that you get the same by just calling the function.


def test_calculate_wilson_outliers(wilson_test_reflection_table):
    """Test the set wilson outliers function."""
    reflection_table = set_wilson_outliers(wilson_test_reflection_table)

    assert list(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling)
    ) == [True, False, True, False]


def test_calculate_prescaling_correction():
    """Test the helper function that applies the lp, dqe and partiality corr."""
    reflection_table = flex.reflection_table()
    reflection_table["lp"] = flex.double([1.0, 0.9, 0.8])
    reflection_table["qe"] = flex.double([0.6, 0.5, 0.4])

    reflection_table = calculate_prescaling_correction(reflection_table)
    assert list(reflection_table["prescaling_correction"]) == [
        1.0 / 0.6,
        0.9 / 0.5,
        0.8 / 0.4,
    ]

    # Test compatibilty for old datasets
    del reflection_table["qe"]
    reflection_table["dqe"] = flex.double([0.6, 0.5, 0.4])
    reflection_table = calculate_prescaling_correction(reflection_table)
    assert list(reflection_table["prescaling_correction"]) == [
        1.0 / 0.6,
        0.9 / 0.5,
        0.8 / 0.4,
    ]


def test_reasons():
    """Test the reasons class, which is basically a dictionary with a nice
    printed output"""
    reasons = Reasons()
    reasons.add_reason("test reason", 100)
    assert reasons.reasons["test reason"] == 100
    print(reasons)
    expected_output = """Reflections passing individual criteria:
criterion: test reason, reflections: 100
"""
    assert reasons.__repr__() == expected_output


def test_calculate_harmonic_tables_from_selections():
    selection = flex.size_t([1, 0, 2, 3, 1])
    coefficients = [flex.double([10, 11, 12, 13]), flex.double([20, 21, 22, 23])]

    arrays, mat = calculate_harmonic_tables_from_selections(
        selection, selection, coefficients
    )
    assert len(arrays) == 2
    assert mat.n_cols == 2
    assert mat.n_rows == 5
    assert list(arrays[0]) == [11, 10, 12, 13, 11]
    assert list(arrays[1]) == [21, 20, 22, 23, 21]
    assert mat[0, 0] == arrays[0][0]
    assert mat[1, 0] == arrays[0][1]
    assert mat[2, 0] == arrays[0][2]
    assert mat[3, 0] == arrays[0][3]
    assert mat[4, 0] == arrays[0][4]
    assert mat[0, 1] == arrays[1][0]
    assert mat[1, 1] == arrays[1][1]
    assert mat[2, 1] == arrays[1][2]
    assert mat[3, 1] == arrays[1][3]
    assert mat[4, 1] == arrays[1][4]


def test_equality_of_two_harmonic_table_methods(dials_regression, run_in_tmpdir):
    from dials_scaling_ext import calc_theta_phi, calc_lookup_index
    from dxtbx.serialize import load
    from dials.util.options import OptionParser
    from libtbx import phil
    from dials.algorithms.scaling.scaling_library import create_scaling_model

    data_dir = os.path.join(dials_regression, "xia2-28")
    pickle_path = os.path.join(data_dir, "20_integrated.pickle")
    sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")

    phil_scope = phil.parse(
        """
      include scope dials.command_line.scale.phil_scope
    """,
        process_includes=True,
    )
    optionparser = OptionParser(phil=phil_scope, check_format=False)
    params, _ = optionparser.parse_args(args=[], quick_parse=True)
    params.model = "physical"
    lmax = 2
    params.parameterisation.lmax = lmax

    reflection_table = flex.reflection_table.from_pickle(pickle_path)
    experiments = load.experiment_list(sweep_path, check_format=False)
    experiments = create_scaling_model(params, experiments, [reflection_table])

    experiment = experiments[0]
    # New method

    reflection_table["phi"] = (
        reflection_table["xyzobs.px.value"].parts()[2]
        * experiment.scan.get_oscillation()[1]
    )
    reflection_table = calc_crystal_frame_vectors(reflection_table, experiment)
    theta_phi_0 = calc_theta_phi(reflection_table["s0c"])  # array of tuples in radians
    theta_phi_1 = calc_theta_phi(reflection_table["s1c"])
    points_per_degree = 2
    s0_lookup_index = calc_lookup_index(theta_phi_0, points_per_degree)
    s1_lookup_index = calc_lookup_index(theta_phi_1, points_per_degree)
    print(list(s0_lookup_index[0:20]))
    print(list(s1_lookup_index[0:20]))
    coefficients_list = create_sph_harm_lookup_table(lmax, points_per_degree)
    experiment.scaling_model.components["absorption"].data = {
        "s0_lookup": s0_lookup_index,
        "s1_lookup": s1_lookup_index,
    }
    experiment.scaling_model.components[
        "absorption"
    ].coefficients_list = coefficients_list
    assert experiment.scaling_model.components["absorption"]._mode == "memory"
    experiment.scaling_model.components["absorption"].update_reflection_data()
    absorption = experiment.scaling_model.components["absorption"]
    harmonic_values_list = absorption.harmonic_values[0]

    experiment.scaling_model.components["absorption"].parameters = flex.double(
        [0.1, -0.1, 0.05, 0.02, 0.01, -0.05, 0.12, -0.035]
    )
    scales, derivatives = experiment.scaling_model.components[
        "absorption"
    ].calculate_scales_and_derivatives()

    # Old method:

    old_data = {"sph_harm_table": create_sph_harm_table(theta_phi_0, theta_phi_1, lmax)}
    experiment.scaling_model.components["absorption"].data = old_data
    assert experiment.scaling_model.components["absorption"]._mode == "speed"
    experiment.scaling_model.components["absorption"].update_reflection_data()
    old_harmonic_values = absorption.harmonic_values[0]
    for i in range(0, 8):
        print(i)
        assert list(harmonic_values_list[i]) == pytest.approx(
            list(old_harmonic_values.col(i).as_dense_vector()), abs=0.01
        )

    experiment.scaling_model.components["absorption"].parameters = flex.double(
        [0.1, -0.1, 0.05, 0.02, 0.01, -0.05, 0.12, -0.035]
    )
    scales_1, derivatives_1 = experiment.scaling_model.components[
        "absorption"
    ].calculate_scales_and_derivatives()

    assert list(scales_1) == pytest.approx(list(scales), abs=0.001)
    assert list(scales_1) != [1.0] * len(scales_1)


def test_calc_lookup_index():
    pi = 3.141
    theta_phi = flex.vec2_double([(0.001, -pi), (pi, pi), (0.001, pi), (pi, -pi)])
    indices = calc_lookup_index(theta_phi, 1)
    assert list(indices) == [0, 64799, 359, 64440]
    indices = calc_lookup_index(theta_phi, 2)
    assert list(indices) == [0, 259199, 719, 258480]
