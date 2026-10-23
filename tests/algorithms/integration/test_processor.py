from __future__ import annotations

import math
import os
from unittest import mock

import psutil
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

import dials.algorithms.integration.processor
from dials.algorithms.profile_model.gaussian_rs import Model
from dials.array_family import flex
from dials_algorithms_integration_integrator_ext import JobList, max_memory_needed


def test_shoebox_memory_is_a_reasonable_guesstimate(dials_data):
    path = dials_data("centroid_test_data") / "experiments.json"

    exlist = ExperimentListFactory.from_json_file(path)[0]
    exlist.profile = Model(
        None,
        n_sigma=3,
        sigma_b=0.024 * math.pi / 180.0,
        sigma_m=0.044 * math.pi / 180.0,
    )

    rlist = flex.reflection_table.from_predictions(exlist)
    rlist["id"] = flex.int(len(rlist), 0)
    rlist["bbox"] = flex.int6(rlist.size(), (0, 1, 0, 1, 0, 1))

    jobs = JobList()
    jobs.add((0, 1), (0, 9), 9, 0)
    # Every single-pixel shoebox ends on frame 0, so all 1996 of them hold
    # data, mask and background (9 bytes per pixel) at the same time, on top
    # of the Shoebox object and one frame lookup index (8 bytes) each
    expected = 1996 * (9 + 8 + _shoebox_object_bytes())
    for flatten in (True, False):
        assumed_memory_usage = list(jobs.shoebox_memory(rlist, flatten))
        assert assumed_memory_usage == [expected]

        max_mem = max_memory_needed(rlist, 0, 9, flatten)
        assert max_mem == expected


def _reflections(bboxes):
    rlist = flex.reflection_table()
    rlist["id"] = flex.int(len(bboxes), 0)
    rlist["flags"] = flex.size_t(len(bboxes), 0)
    rlist["bbox"] = flex.int6(bboxes)
    return rlist


def _shoebox_object_bytes():
    """The fixed cost of an unallocated Shoebox object, which is derived
    from sizeof() so depends on the platform"""
    # A single pixel on a single frame costs 9 bytes of pixel data, plus one
    # frame lookup index
    return max_memory_needed(_reflections([(0, 1, 0, 1, 0, 1)]), 0, 1, False) - 17


def test_shoebox_object_bytes_is_plausible():
    # Measured as 448 bytes per unallocated shoebox on 64-bit Linux
    assert 300 < _shoebox_object_bytes() < 600


@pytest.mark.parametrize(
    "flatten,expected_pixel_memory",
    [
        # Peak is on the last frame of A: data and mask (5 bytes per pixel) for
        # all 4 frames of A, plus background (4 bytes per pixel). B has been freed.
        (False, 5 * 400 + 4 * 400),
        # Flattened shoeboxes have a single frame of 100 pixels for A
        (True, 5 * 100 + 4 * 100),
    ],
)
def test_shoebox_memory_counts_background_only_on_final_frame(
    flatten, expected_pixel_memory
):
    rlist = _reflections(
        [
            (0, 10, 0, 10, 0, 4),  # A: 10x10 pixels on frames 0-3
            (0, 5, 0, 5, 1, 2),  # B: 5x5 pixels on frame 1
        ]
    )
    # Both Shoebox objects and the frame lookup (4 + 1 entries) exist for the
    # whole job
    expected = expected_pixel_memory + 2 * _shoebox_object_bytes() + 5 * 8

    jobs = JobList()
    jobs.add((0, 1), (0, 4), 4, 0)
    assert list(jobs.shoebox_memory(rlist, flatten)) == [expected]
    assert max_memory_needed(rlist, 0, 4, flatten) == expected


def test_shoebox_memory_includes_bytes_per_reflection():
    rlist = _reflections([(0, 10, 0, 10, 0, 4), (0, 5, 0, 5, 1, 2)])
    jobs = JobList()
    jobs.add((0, 1), (0, 4), 4, 0)
    for flatten in (True, False):
        base = max_memory_needed(rlist, 0, 4, flatten)
        assert max_memory_needed(rlist, 0, 4, flatten, 100) == base + 2 * 100
        assert list(jobs.shoebox_memory(rlist, flatten, 100)) == [base + 2 * 100]


@mock.patch("dials.algorithms.integration.processor.flex.max")
@mock.patch("dials.algorithms.integration.processor.MEMORY_LIMIT", 1_000_000_000)
def test_runtime_error_raised_when_not_enough_memory(mock_flex_max):
    # We are taking into account current memory usage, which will typically be ~200MB.
    # So set the shoebox memory limit to be 500MB, which should fail the memory test.
    mock_flex_max.return_value = 500_000_000

    phil_mock = mock.Mock()
    phil_mock.mp.method = "multiprocessing"
    phil_mock.mp.nproc = 4
    phil_mock.block.max_memory_usage = 0.5

    reflections = {"bbox": flex.int6(1000, (0, 1, 0, 1, 0, 1))}
    # No detectors, so no memory is required for images
    experiments = mock.Mock()
    experiments.detectors.return_value = []
    manager = dials.algorithms.integration.processor._Manager(
        experiments, reflections, phil_mock
    )
    manager.jobs = mock.Mock(autospec=JobList)

    with pytest.raises(MemoryError) as exc_info:
        manager.compute_processors()
    assert "Not enough memory to run integration jobs." in exc_info.value.args[0]
    mock_flex_max.assert_called_once_with(manager.jobs.shoebox_memory.return_value)

    # Reduce memory usage of the shoeboxes by the current memory usage
    # (with an extra 1MB for safety), should then pass.
    current_memory_usage = psutil.Process(os.getpid()).memory_info().rss
    mock_flex_max.return_value = 499_000_000 - current_memory_usage
    manager.compute_processors()
    mock_flex_max.assert_called_with(manager.jobs.shoebox_memory.return_value)
