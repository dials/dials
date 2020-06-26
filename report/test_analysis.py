"""Tests for dials.report.analysis module"""
from __future__ import absolute_import, division, print_function
import pytest
from mock import mock
from cctbx import miller
from dxtbx.model import Crystal
from dxtbx.serialize import load
from dials.array_family import flex
from dials.report.analysis import (
    scales_vs_batch,
    i_sig_i_vs_batch,
    rmerge_vs_batch,
    batch_dependent_properties,
    reflection_tables_to_batch_dependent_properties,
    combined_table_to_batch_dependent_properties,
    make_xia2_style_statistics_summary,
)
from dials.algorithms.scaling.scaling_library import (
    merging_stats_from_scaled_array,
    scaled_data_as_miller_array,
)


@pytest.fixture
def example_crystal():
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 1.0],
        "space_group_hall_symbol": "-P 2yb",
    }
    crystal = Crystal.from_dict(exp_dict)
    return crystal


@pytest.fixture
def example_miller_set(example_crystal):
    """Generate an example miller set."""
    ms = miller.set(
        crystal_symmetry=example_crystal.get_crystal_symmetry(),
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


def test_reflections_to_batch_properties(
    data_array, example_miller_set, example_crystal
):
    """Test the helper functions that provide the batch properties from reflection
    tables and experiments."""
    # first make a reflection table.
    reflections = flex.reflection_table()
    reflections["intensity.scale.value"] = data_array.data() * flex.double(9, 2.0)
    reflections["inverse_scale_factor"] = flex.double(9, 2.0)
    reflections["intensity.scale.variance"] = flex.double(9, 4.0) * flex.double(9, 4.0)
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [(0, 0, 0.1)] * 3 + [(0, 0, 1.1)] * 3 + [(0, 0, 2.1)] * 2 + [(0, 0, 3.1)]
    )
    reflections["miller_index"] = example_miller_set.indices()
    reflections["id"] = flex.int(9, 1)
    reflections.set_flags(flex.bool(9, True), reflections.flags.integrated)

    experiments = [mock.Mock()]
    experiments[0].scan.get_image_range.return_value = [1, 10]
    experiments[0].crystal = example_crystal

    (
        bins,
        rmerge,
        isigi,
        scalesvsbatch,
        batch_data,
    ) = reflection_tables_to_batch_dependent_properties(  # pylint: disable=unbalanced-tuple-unpacking
        [reflections], experiments
    )

    assert bins == expected_results["bins"]
    assert rmerge == pytest.approx(expected_results["rmergevb"], 1e-4)
    assert isigi == pytest.approx(expected_results["isigivb"], 1e-4)
    assert scalesvsbatch == pytest.approx([2.0] * 4, 1e-4)
    assert batch_data == [{"range": (1, 10), "id": 0}]

    # now try a two experiment dataset in a combined table.
    import copy

    reflections_2 = copy.deepcopy(reflections)
    reflections_2["id"] = flex.int(9, 2)
    reflections.extend(reflections_2)
    experiments = [mock.Mock(), mock.Mock()]
    experiments[0].scan.get_image_range.return_value = [1, 10]
    experiments[0].crystal = example_crystal
    experiments[1].scan.get_image_range.return_value = [1, 10]
    experiments[1].crystal = example_crystal

    (
        bins,
        rmerge,
        isigi,
        scalesvsbatch,
        batch_data,
    ) = combined_table_to_batch_dependent_properties(  # pylint: disable=unbalanced-tuple-unpacking
        reflections, experiments
    )

    assert bins == [1, 2, 3, 4, 101, 102, 103, 104]
    assert rmerge == pytest.approx(expected_results["rmergevb"] * 2, 1e-4)
    assert isigi == pytest.approx(expected_results["isigivb"] * 2, 1e-4)
    assert scalesvsbatch == pytest.approx([2.0] * 8, 1e-4)
    assert batch_data == [{"range": (1, 10), "id": 0}, {"range": (101, 110), "id": 1}]


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


def test_make_xia2_style_statistics_summary(dials_data):

    location = dials_data("l_cysteine_4_sweeps_scaled")
    expts = load.experiment_list(location.join("scaled_20_25.expt"), check_format=False)
    refls = flex.reflection_table.from_file(location.join("scaled_20_25.refl"))
    # Get a miller array of real data and calculate an iotbx.merging_statistics
    ma = scaled_data_as_miller_array([refls], expts)
    arr, anom = merging_stats_from_scaled_array(ma)

    # Test that something is returned in each case
    ### Case of overall statistics summary
    out = make_xia2_style_statistics_summary(arr, anom)
    assert out
    assert all(a in out for a in ("Overall", "Low", "High"))
    assert "Suggested" not in out
    ### Case of overall and suggested statistics summary (with anom)
    out = make_xia2_style_statistics_summary(arr, anom, arr, anom)
    assert out
    assert all(a in out for a in ("Overall", "Suggested", "Low", "High"))
    ### Case of no anomalous, but with suggested as well as overall.
    out = make_xia2_style_statistics_summary(arr, selected_statistics=arr)
    assert out
    assert all(a in out for a in ("Overall", "Suggested", "Low", "High"))
