"""Tests for dials.report.analysis module"""
import pytest
from cctbx import miller
from cctbx.array_family import flex
from dxtbx.model import Crystal
from dials.report.analysis import (
    scales_vs_batch,
    i_sig_i_vs_batch,
    rmerge_vs_batch,
    batch_dependent_properties,
)


@pytest.fixture
def example_miller_set():
    """Generate an example miller set."""
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 1.0],
        "space_group_hall_symbol": "-P 2yb",
    }
    crystal = Crystal.from_dict(exp_dict)
    ms = miller.set(
        crystal_symmetry=crystal.get_crystal_symmetry(),
        indices=flex.miller_index([(1, 1, 1)] * 8 + [(2, 2, 2)]),
        anomalous_flag=False,
    )
    return ms


@pytest.fixture
def batch_array(example_miller_set):
    """Generate an example batch array."""
    batches_ = flex.double([1, 1, 1, 2, 2, 2, 3, 3, 4])
    batches = miller.array(example_miller_set, data=batches_)
    return batches


@pytest.fixture
def data_array(example_miller_set):
    """Generate an example data array."""
    data_ = flex.double([1.0 + (i * 0.1) for i in range(0, 9)])
    data = miller.array(example_miller_set, data=data_)
    return data


expected_results = {
    "bins": [1, 2, 3, 4],
    "svb": [1.1, 1.4, 1.65, 1.8],
    "isigivb": [1.1 / 2.0, 1.4 / 2.0, 1.65 / 2.0, 1.8 / 2.0],
    "rmergevb": [0.22727, 0.059524, 0.18182, 0.0],
}


def test_scales_vs_batch(batch_array, data_array):
    """Test the scales_vs_batch function."""
    bins, svb = scales_vs_batch(data_array, batch_array)
    assert bins == expected_results["bins"]
    assert svb == pytest.approx(expected_results["svb"], 1e-6)


def test_IsigI_vs_batch(batch_array, data_array):
    """Test the IsigI_vs_batch function."""
    with pytest.raises(AssertionError):
        bins, isigivb = i_sig_i_vs_batch(data_array, batch_array)

    Is = data_array.customized_copy()
    Is.set_sigmas(flex.double(9, 2.0))
    bins, isigivb = i_sig_i_vs_batch(Is, batch_array)
    assert bins == expected_results["bins"]
    assert isigivb == pytest.approx(expected_results["isigivb"], 1e-6)


def test_Rmerge_vs_batch(batch_array, data_array):
    """Test the Rmerge_vs_batch function."""
    bins, rmergevsb = rmerge_vs_batch(data_array, batch_array)
    assert bins == expected_results["bins"]
    assert rmergevsb == pytest.approx(expected_results["rmergevb"], 1e-4)


def test_batch_dependent_properties(batch_array, data_array):
    """Test the interface function that manages the calculations."""
    Is = data_array.customized_copy()
    Is.set_sigmas(flex.double(9, 2.0))

    bins, rmerge, isigi, scalesvsbatch = batch_dependent_properties(
        batch_array, Is, scales=data_array
    )

    assert bins == expected_results["bins"]
    assert rmerge == pytest.approx(expected_results["rmergevb"], 1e-4)
    assert isigi == pytest.approx(expected_results["isigivb"], 1e-4)
    assert scalesvsbatch == pytest.approx(expected_results["svb"], 1e-4)

    # test for no scales given
    bins, rmerge, isigi, scalesvsbatch = batch_dependent_properties(batch_array, Is)

    assert bins == expected_results["bins"]
    assert rmerge == pytest.approx(expected_results["rmergevb"], 1e-4)
    assert isigi == pytest.approx(expected_results["isigivb"], 1e-4)
    assert scalesvsbatch is None

    # test for bad input - batches should reduce for all data
    Is.set_sigmas(flex.double([2.0] * 8 + [-1.0] * 1))
    bins, rmerge, isigi, scalesvsbatch = batch_dependent_properties(
        batch_array, Is, scales=data_array
    )

    assert bins == expected_results["bins"][0:3]
    assert rmerge == pytest.approx(expected_results["rmergevb"][0:3], 1e-4)
    assert isigi == pytest.approx(expected_results["isigivb"][0:3], 1e-4)
    assert scalesvsbatch == pytest.approx(expected_results["svb"][0:3], 1e-4)

    bins, rmerge, isigi, scalesvsbatch = batch_dependent_properties(batch_array, Is)

    assert bins == expected_results["bins"][0:3]
    assert rmerge == pytest.approx(expected_results["rmergevb"][0:3], 1e-4)
    assert isigi == pytest.approx(expected_results["isigivb"][0:3], 1e-4)
    assert scalesvsbatch is None

    # test for mismatched array sizes
    with pytest.raises(AssertionError):
        _ = batch_dependent_properties(batch_array, Is, data_array[0:-1])
    with pytest.raises(AssertionError):
        _ = batch_dependent_properties(batch_array[0:-1], Is)
