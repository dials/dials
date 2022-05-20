from __future__ import annotations

import math
import os
from unittest import mock

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

import dials.algorithms.integration.processor
from dials.algorithms.integration.processor import assess_available_memory
from dials.algorithms.profile_model.gaussian_rs import Model
from dials.array_family import flex
from dials_algorithms_integration_integrator_ext import JobList, max_memory_needed


def test_shoebox_memory_is_a_reasonable_guesstimate(dials_data):
    path = dials_data("centroid_test_data", pathlib=True) / "experiments.json"

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
@mock.patch("dials.algorithms.integration.processor.psutil.virtual_memory")
@mock.patch("dials.algorithms.integration.processor.psutil.swap_memory")
def test_runtime_error_raised_when_not_enough_memory(
    mock_psutil_swap, mock_psutil_vm, mock_flex_max
):
    mock_flex_max.return_value = 750001
    mock_psutil_vm.return_value.available = 1000000
    mock_psutil_swap.return_value.free = 0

    phil_mock = mock.Mock()
    phil_mock.mp.method = "multiprocessing"
    phil_mock.mp.nproc = 4
    phil_mock.block.max_memory_usage = 0.75

    reflections = {"bbox": flex.int6(1000, (0, 1, 0, 1, 0, 1))}
    manager = dials.algorithms.integration.processor._Manager(
        None, reflections, phil_mock
    )
    manager.jobs = mock.Mock(autospec=JobList)

    with pytest.raises(MemoryError) as exc_info:
        manager.compute_processors()
    assert "Not enough memory to run integration jobs." in exc_info.value.args[0]
    mock_flex_max.assert_called_once_with(manager.jobs.shoebox_memory.return_value)

    # Reduce memory usage by 1 byte, should then pass
    mock_flex_max.return_value = 750000
    manager.compute_processors()
    mock_flex_max.assert_called_with(manager.jobs.shoebox_memory.return_value)


def test_assess_available_memory_condor_job_ad(mocker, monkeypatch, tmp_path):
    job_ad = tmp_path / ".job.ad"
    monkeypatch.setenv("_CONDOR_JOB_AD", os.fspath(job_ad))

    # Test that we handle gracefully absence of file
    params = mocker.Mock()
    params.block.max_memory_usage = 0.9
    available_memory, _, _ = assess_available_memory(params)
    assert available_memory and available_memory != 123

    # Now write something sensible to the file
    job_ad.write_text(
        """\
MemoryProvisioned = 123
"""
    )
    available_memory, _, _ = assess_available_memory(params)
    assert available_memory == params.block.max_memory_usage * 123 * 1024**2

    # Now check it is robust against parsing failures
    job_ad.write_text(
        """\
MemoryProvisioned =
"""
    )
    available_memory, _, _ = assess_available_memory(params)
    assert available_memory and available_memory != 123
