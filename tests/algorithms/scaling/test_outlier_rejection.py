"""
Tests for outlier rejection.
"""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from cctbx import uctbx
from cctbx.sgtbx import space_group

from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.outlier_rejection import (
    NormDevOutlierRejection,
    SimpleNormDevOutlierRejection,
    TargetedOutlierRejection,
    determine_Esq_outlier_index_arrays,
    determine_outlier_index_arrays,
    limit_outlier_weights,
    reject_outliers,
)
from dials.array_family import flex


@pytest.fixture(scope="module")
def test_sg():
    """Create a space group object."""
    return space_group("C 2y")


@pytest.fixture
def mock_exp_with_sg(test_sg):
    """Create a mock experiment with a space group"""
    exp = Mock()
    exp.crystal.get_space_group.return_value = test_sg
    exp.crystal.get_unit_cell.return_value = uctbx.unit_cell((5, 5, 5, 90, 90, 90))
    return exp


@pytest.fixture
def generated_Ih_table(test_sg):
    """Generate an Ih_table"""
    rt = generate_outlier_table()
    Ih_table = IhTable([rt], test_sg, nblocks=1)
    return Ih_table


@pytest.fixture
def outlier_target_table(test_sg):
    """Generate an Ih_table for targeted outlier rejection"""
    target = generate_target_table()
    target_Ih = IhTable([target], test_sg, nblocks=1)
    return target_Ih


def generate_target_table():
    """Reflection table defining target reflections"""
    target = flex.reflection_table()
    target["intensity"] = flex.double([500])
    target["variance"] = flex.double([1.0])
    target["inverse_scale_factor"] = flex.double([1.0])
    target["miller_index"] = flex.miller_index([(0, 0, 2)])
    target.set_flags(flex.bool([False]), target.flags.excluded_for_scaling)
    target.set_flags(flex.bool([False]), target.flags.user_excluded_in_scaling)
    return target


def generate_outlier_table():
    """Generate a reflection table for outlier testing."""
    rt = flex.reflection_table()
    rt["intensity"] = flex.double(
        [1.0, 1.0, 1.0, 1.0, 20.0, 1.0, 1.0, 20.0, 400.0, 500.0, 10.0, 10.0]
    )
    rt["variance"] = flex.double(12, 1.0)
    rt["inverse_scale_factor"] = flex.double(12, 1.0)
    rt["miller_index"] = flex.miller_index(
        [
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 23),
            (0, 0, 3),
        ]
    )
    rt.set_flags(
        flex.bool(
            [
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ]
        ),
        rt.flags.excluded_for_scaling,
    )
    rt.set_flags(flex.bool(12, False), rt.flags.user_excluded_in_scaling)
    return rt


def generate_outlier_table_2():
    """Generate a reflection table for outlier testing."""
    rt = flex.reflection_table()
    rt["intensity"] = flex.double(
        [1.0, 1.0, 1.0, 20.0, 1.0, 1.0, 20.0, 400.0, 500.0, 10.0, 10.0]
    )
    rt["variance"] = flex.double(11, 1.0)
    rt["miller_index"] = flex.miller_index(
        [
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 23),
            (0, 0, 3),
        ]
    )
    rt.set_flags(
        flex.bool(
            [
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ]
        ),
        rt.flags.excluded_for_scaling,
    )
    return rt


expected_standard_output = [
    False,
    False,
    False,
    False,
    True,
    False,
    False,
    True,
    True,
    True,
    False,
    False,
]

expected_simple_output = [
    False,
    False,
    False,
    False,
    True,
    True,
    True,
    True,
    True,
    True,
    False,
    False,
]

expected_target_output = [
    False,
    False,
    False,
    False,
    False,
    True,
    True,
    True,
    True,
    False,
    False,
    False,
]


def test_outlier_rejection_with_small_outliers():

    rt = flex.reflection_table()
    rt["intensity"] = flex.double(
        [3560.84231, 3433.66407, 3830.64235, 0.20552, 3786.59537]
        + [4009.98652, 0.00000, 3578.91470, 3549.19151, 3379.58616]
        + [3686.38610, 3913.42869, 0.00000, 3608.84869, 3681.11110]
    )
    rt["variance"] = flex.double(
        [10163.98104, 9577.90389, 9702.84868, 3.77427, 8244.70685]
        + [9142.38221, 1.51118, 9634.53782, 9870.73103, 9078.23488]
        + [8977.26984, 8712.91360, 1.78802, 7473.26521, 10075.49862]
    )
    rt["inverse_scale_factor"] = flex.double(rt.size(), 1.0)
    rt["miller_index"] = flex.miller_index([(0, 0, 1)] * rt.size())
    expected_outliers = [3, 6, 12]

    OutlierRej = NormDevOutlierRejection(IhTable([rt], space_group("P 1")), zmax=6.0)
    OutlierRej.run()
    outliers = OutlierRej.final_outlier_arrays
    assert len(outliers) == 1
    assert set(outliers[0]) == set(expected_outliers)


def test_limit_outlier_weights():

    rt = flex.reflection_table()
    rt["intensity"] = flex.double([100.0, 101.0, 109.0, 105.0, 1.0])
    rt["variance"] = flex.double([100.0, 101.0, 109.0, 105.0, 1.0])

    rt["inverse_scale_factor"] = flex.double(rt.size(), 1.0)
    rt["miller_index"] = flex.miller_index([(0, 0, 1)] * rt.size())
    rt2 = flex.reflection_table()
    rt2["intensity"] = flex.double([100.0, 101.0, 102.0, 105.0, 1.0])
    rt2["variance"] = flex.double([100.0, 101.0, 102.0, 105.0, 1.0])
    rt2["inverse_scale_factor"] = flex.double(rt.size(), 1.0)
    rt2["miller_index"] = flex.miller_index([(0, 0, 1)] * rt.size())

    table = IhTable([rt, rt2], space_group("P 1"))
    import copy

    new_weights = limit_outlier_weights(
        copy.deepcopy(table.Ih_table_blocks[0].weights),
        table.Ih_table_blocks[0].h_index_matrix,
    )
    assert all(i <= 0.1 for i in new_weights)


def test_multi_dataset_outlier_rejection(test_sg):
    """Test outlier rejection with two datasets."""
    rt1 = flex.reflection_table()
    rt1["intensity"] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 20.0, 400.0, 10.0])
    rt1["variance"] = flex.double(9, 1.0)
    rt1["inverse_scale_factor"] = flex.double(9, 1.0)
    rt1["miller_index"] = flex.miller_index(
        [
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 2),
            (0, 0, 3),
        ]
    )
    rt1.set_flags(
        flex.bool([True, False, False, False, False, False, False, False, False]),
        rt1.flags.excluded_for_scaling,
    )
    rt1.set_flags(flex.bool(9, False), rt1.flags.user_excluded_in_scaling)
    rt2 = flex.reflection_table()
    rt2["intensity"] = flex.double([10.0, 20.0, 500.0])
    rt2["variance"] = flex.double(3, 1.0)
    rt2["inverse_scale_factor"] = flex.double(3, 1.0)
    rt2["miller_index"] = flex.miller_index([(0, 0, 23), (0, 0, 1), (0, 0, 2)])
    rt2.set_flags(flex.bool([False, False, False]), rt1.flags.excluded_for_scaling)
    rt2.set_flags(flex.bool(3, False), rt2.flags.user_excluded_in_scaling)
    Ih_table = IhTable([rt1, rt2], test_sg, nblocks=1)
    zmax = 6.0
    outlier_rej = SimpleNormDevOutlierRejection(Ih_table, zmax)
    outlier_rej.run()
    outliers = outlier_rej.final_outlier_arrays
    assert len(outliers) == 2
    assert list(outliers[0]) == [4, 5, 6, 7]
    assert list(outliers[1]) == [1, 2]

    outlier_rej = NormDevOutlierRejection(Ih_table, zmax)
    outlier_rej.run()
    outliers = outlier_rej.final_outlier_arrays
    assert len(outliers) == 2
    assert list(outliers[0]) == [7, 6]
    assert list(outliers[1]) == [1, 2]


def test_standard_outlier_rejection(generated_Ih_table):
    """Test the standard outlier rejection algorithm."""
    zmax = 6.0
    OutlierRej = NormDevOutlierRejection(generated_Ih_table, zmax)
    OutlierRej.run()
    outliers = OutlierRej.final_outlier_arrays
    assert len(outliers) == 1
    assert list(outliers[0]) == [4, 9, 8, 7]


def test_targeted_outlier_rejection(generated_Ih_table, outlier_target_table):
    """Test the targeted outlier rejection algorithm - only reflections
    that exist in both the target and the reflecton table should be tested
    based on their normalised deviation."""
    zmax = 6.0
    outlier_rej = TargetedOutlierRejection(
        generated_Ih_table, zmax, outlier_target_table
    )
    outlier_rej.run()
    outliers = outlier_rej.final_outlier_arrays
    assert len(outliers) == 1
    assert list(outliers[0]) == [5, 6, 7, 8]


def test_simple_outlier_rejection(generated_Ih_table):
    """Test the simple outlier rejection algorithm."""
    zmax = 6.0
    outlier_rej = SimpleNormDevOutlierRejection(generated_Ih_table, zmax)
    outlier_rej.run()
    outliers = outlier_rej.final_outlier_arrays
    assert len(outliers) == 1
    assert list(outliers[0]) == [4, 5, 6, 7, 8, 9]


def test_reject_outliers(mock_exp_with_sg):
    """Test the reject outliers function"""
    refls = generate_outlier_table_2()
    exp = mock_exp_with_sg
    refls = reject_outliers(refls, exp)
    assert (
        list(refls.get_flags(refls.flags.outlier_in_scaling))
        == expected_standard_output[1:]
    )

    # Try the case for two tables joined together
    refls = generate_outlier_table_2()
    refls2 = generate_outlier_table_2()
    refls.extend(refls2)
    refls = reject_outliers(refls, exp)
    assert (
        list(refls.get_flags(refls.flags.outlier_in_scaling))
        == expected_standard_output[1:] * 2
    )

    # Check that any existing outlier flags are overwritten
    refls = generate_outlier_table_2()
    refls.set_flags(flex.bool([True] + [False] * 10), refls.flags.outlier_in_scaling)
    refls = reject_outliers(refls, exp)
    assert (
        list(refls.get_flags(refls.flags.outlier_in_scaling))
        == expected_standard_output[1:]
    )

    # Test for suitable error raising if intensity column not present
    del refls["intensity"]
    with pytest.raises(AssertionError):
        refls = reject_outliers(refls, exp)


def test_determine_outlier_index_arrays(generated_Ih_table, outlier_target_table):
    """Test the helper function."""

    outliers = determine_outlier_index_arrays(generated_Ih_table, "standard")
    assert list(outliers[0]) == [4, 9, 8, 7]

    outliers = determine_outlier_index_arrays(generated_Ih_table, "simple")
    assert list(outliers[0]) == [4, 5, 6, 7, 8, 9]

    outliers = determine_outlier_index_arrays(
        generated_Ih_table, "target", target=outlier_target_table
    )
    assert list(outliers[0]) == [5, 6, 7, 8]

    outliers = determine_outlier_index_arrays(generated_Ih_table, method=None)
    assert len(outliers) == 1
    assert not outliers[0]
    with pytest.raises(ValueError):
        _ = determine_outlier_index_arrays(generated_Ih_table, "badchoice")[0]


def test_determine_Esq_outlier_index_arrays(
    generated_Ih_table, mock_exp_with_sg, test_sg
):
    # Set the emax lower to check that two reflections are identified as outliers
    outliers = determine_Esq_outlier_index_arrays(
        generated_Ih_table, mock_exp_with_sg, emax=1.5
    )
    assert list(outliers[0]) == [8, 9]
    # now split the dataset into two, to check the output is correctly formed
    rt = generate_outlier_table()
    rt1 = rt[0:9]
    rt2 = rt[9:]
    Ih_table = IhTable([rt1, rt2], test_sg)
    outliers = determine_Esq_outlier_index_arrays(Ih_table, mock_exp_with_sg, emax=1.5)
    assert list(outliers[0]) == [8]
    assert list(outliers[1]) == [0]
    assert len(outliers) == 2
