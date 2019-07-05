"""
Test for the scaling weighting module.
"""
from __future__ import absolute_import, division, print_function
import pytest
from mock import Mock
from dials.array_family import flex
from dials.algorithms.scaling.weighting import (
    get_weighting_scheme,
    WeightingBase,
    UnityWeights,
    GemanMcClureWeights,
    HuberWeights,
    CauchyWeights,
)


@pytest.fixture
def mock_Ih_table():
    """Mock of an Ih table for testing the weighting schemes."""
    Ih_table = Mock()
    Ih_table.variances = flex.double([1.0, 1.0, 2.0, 2.0])
    Ih_table.intensities = flex.double([1.0, 1.0, 2.0, 4.5])
    Ih_table.inverse_scale_factors = flex.double(4, 1.0)
    Ih_table.Ih_values = flex.double(4, 1.5)
    Ih_table.size = 4
    return Ih_table


def test_weighting_schemes(mock_Ih_table):
    """Test for the different weighting scheme classes and helper function."""

    weights = get_weighting_scheme(mock_Ih_table, "invvar")
    assert isinstance(weights, WeightingBase)
    assert weights.weighting_scheme == "invvar"
    weights.calculate_initial_weights()
    assert list(mock_Ih_table.weights) == [1.0, 1.0, 0.5, 0.5]

    weights = get_weighting_scheme(mock_Ih_table, "unity")
    assert isinstance(weights, UnityWeights)
    assert weights.weighting_scheme == "unity"
    weights.calculate_initial_weights()
    assert list(mock_Ih_table.weights) == [1.0, 1.0, 1.0, 1.0]

    weights = get_weighting_scheme(mock_Ih_table, "GM")
    assert isinstance(weights, GemanMcClureWeights)
    assert weights.weighting_scheme == "GM"
    weights.calculate_initial_weights()
    assert list(mock_Ih_table.weights) == [1.0, 1.0, 0.5, 0.5]
    weights.apply_iterative_weights()
    # weights are 1/(1+e^2)^2, e=1/3 for all in the test example, except
    # last which is 2.
    assert list(mock_Ih_table.weights) == pytest.approx(
        [81.0 / 100.0] * 3 + [1.0 / 25.0]
    )

    weights = get_weighting_scheme(mock_Ih_table, "huber")
    assert isinstance(weights, HuberWeights)
    assert weights.weighting_scheme == "huber"
    weights.calculate_initial_weights()
    assert list(mock_Ih_table.weights) == [1.0, 1.0, 0.5, 0.5]
    weights.apply_iterative_weights()
    assert list(mock_Ih_table.weights) == pytest.approx([1.0] * 3 + [weights.c / 2.0])

    weights = get_weighting_scheme(mock_Ih_table, "cauchy")
    assert isinstance(weights, CauchyWeights)
    assert weights.weighting_scheme == "cauchy"
    weights.calculate_initial_weights()
    assert list(mock_Ih_table.weights) == [1.0, 1.0, 0.5, 0.5]
    weights.apply_iterative_weights()
    # weights are 1/(1+(e/c)^2), e=1/3 for all in the test example, except
    # last which is 2.
    assert list(mock_Ih_table.weights) == pytest.approx(
        [1.0 / (1.0 + (1.0 / (3.0 * weights.c)) ** 2)] * 3
        + [1.0 / (1.0 + (2.0 / weights.c) ** 2)]
    )

    with pytest.raises(ValueError):
        weights = get_weighting_scheme(mock_Ih_table, "badchoice")
