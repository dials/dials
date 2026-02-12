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
    for flatten in (True, False):
        assumed_memory_usage = list(jobs.shoebox_memory(rlist, flatten))
        assert len(assumed_memory_usage) == 1
        assert assumed_memory_usage[0] == pytest.approx(23952, abs=3000)

        max_mem = max_memory_needed(rlist, 0, 9, flatten)
        assert max_mem == pytest.approx(23952, abs=3000)


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
    manager = dials.algorithms.integration.processor._Manager(
        None, reflections, phil_mock
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
