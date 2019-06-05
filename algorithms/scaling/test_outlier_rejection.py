"""
Tests for outlier rejection.
"""
from __future__ import absolute_import, division, print_function
import pytest
from mock import Mock
from cctbx.sgtbx import space_group
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.outlier_rejection import (
    reject_outliers,
    NormDevOutlierRejection,
    SimpleNormDevOutlierRejection,
    determine_outlier_index_arrays,
    TargetedOutlierRejection,
)


@pytest.fixture(scope="module")
def test_sg():
    """Create a space group object."""
    return space_group("C 2y")


@pytest.fixture
def mock_exp_with_sg(test_sg):
    """Create a mock experiment with a space group"""
    exp = Mock()
    exp.crystal.get_space_group.return_value = test_sg
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
    outliers = SimpleNormDevOutlierRejection(Ih_table, zmax).final_outlier_arrays
    assert len(outliers) == 2
    assert list(outliers[0]) == [4, 5, 6, 7]
    assert list(outliers[1]) == [1, 2]

    outliers = NormDevOutlierRejection(Ih_table, zmax).final_outlier_arrays
    assert len(outliers) == 2
    assert list(outliers[0]) == [7, 6]
    assert list(outliers[1]) == [1, 2]


def test_standard_outlier_rejection(generated_Ih_table):
    """Test the standard outlier rejection algorithm."""
    zmax = 6.0
    OutlierRej = NormDevOutlierRejection(generated_Ih_table, zmax)
    outliers = OutlierRej.final_outlier_arrays
    assert len(outliers) == 1
    assert list(outliers[0]) == [4, 9, 8, 7]


def test_targeted_outlier_rejection(generated_Ih_table, outlier_target_table):
    """Test the targeted outlier rejection algorithm - only reflections
    that exist in both the target and the reflecton table should be tested
    based on their normalised deviation."""
    zmax = 6.0
    outliers = TargetedOutlierRejection(
        generated_Ih_table, zmax, outlier_target_table
    ).final_outlier_arrays
    assert len(outliers) == 1
    assert list(outliers[0]) == [5, 6, 7, 8]


def test_simple_outlier_rejection(generated_Ih_table):
    """Test the simple outlier rejection algorithm."""
    zmax = 6.0
    outliers = SimpleNormDevOutlierRejection(
        generated_Ih_table, zmax
    ).final_outlier_arrays
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
