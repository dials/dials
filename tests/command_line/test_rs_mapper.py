from __future__ import annotations

import shutil
import subprocess
from os.path import join

import numpy as np
import pytest

from dxtbx import flumpy
from dxtbx.model.experiment_list import ExperimentListFactory
from iotbx import ccp4_map
from scitbx.array_family import flex

from dials.algorithms.rs_mapper import normalize_voxels, process_tof_experiment_list


def test_rs_mapper(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.rs_mapper"),
            dials_data("centroid_test_data") / "imported_experiments.json",
            'map_file="junk.ccp4"',
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "junk.ccp4").is_file()

    # load results
    m = ccp4_map.map_reader(file_name=str(tmp_path / "junk.ccp4"))
    assert len(m.data) == 7189057
    assert m.header_min == 0.0
    assert flex.min(m.data) == 0.0

    assert m.header_max == 2052.75
    assert flex.max(m.data) == 2052.75

    assert m.header_mean == pytest.approx(0.018924040719866753, abs=1e-6)
    assert flex.mean(m.data) == pytest.approx(0.01892407052218914, abs=1e-6)


def test_multi_panel(dials_data, tmp_path):
    data_dir = dials_data("image_examples")
    image = data_dir / "DLS_I23_germ_13KeV_0001.cbf"

    result = subprocess.run(
        [
            shutil.which("dials.rs_mapper"),
            image,
            'map_file="junk.ccp4"',
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "junk.ccp4").is_file()

    # load results
    m = ccp4_map.map_reader(file_name=str(tmp_path / "junk.ccp4"))
    assert len(m.data) == 7189057
    assert m.header_min == 0.0
    assert flex.min(m.data) == 0.0

    assert m.header_max == 31342.25
    assert flex.max(m.data) == 31342.25

    assert m.header_mean == pytest.approx(0.05911629647016525, abs=1e-6)
    assert flex.mean(m.data) == pytest.approx(0.05911629647016525, abs=1e-6)


def test_masked(dials_data, tmp_path):
    subprocess.run(
        [
            shutil.which("dials.import"),
            dials_data("image_examples") / "dectris_eiger_master.h5",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    result = subprocess.run(
        [shutil.which("dials.rs_mapper"), "imported.expt", "map_file=masked.ccp4"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "masked.ccp4").is_file()

    # Also check the ignore case (regressions here may indicate changes in mask handling)
    result = subprocess.run(
        [
            shutil.which("dials.rs_mapper"),
            "imported.expt",
            "map_file=unmasked.ccp4",
            "ignore_mask=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "unmasked.ccp4").is_file()

    # Load results
    masked = ccp4_map.map_reader(file_name=str(tmp_path / "masked.ccp4"))
    unmasked = ccp4_map.map_reader(file_name=str(tmp_path / "unmasked.ccp4"))

    assert masked.header_max < unmasked.header_max
    assert masked.header_max == pytest.approx(289.11111)
    assert unmasked.header_max == pytest.approx(65535.0)


def test_process_tof_experiment_list(dials_data):
    grid_size = 192
    grid = flex.double(flex.grid(grid_size, grid_size, grid_size), 0)
    counts = flex.int(flex.grid(grid_size, grid_size, grid_size), 0)
    max_resolution = 2
    nproc = 1
    image_file = join(dials_data("isis_sxd_example_data"), "sxd_nacl_run.nxs")
    experiments = ExperimentListFactory.from_filenames([image_file])

    record_all_counts = False
    process_tof_experiment_list(
        experiments, max_resolution, grid, counts, nproc, record_all_counts
    )
    normalize_voxels(grid, counts)
    arr = flumpy.to_numpy(grid)
    arr_counts = flumpy.to_numpy(counts)
    assert np.sum(arr) == pytest.approx(91542.09445706864)
    assert np.max(arr) == pytest.approx(1064.0)
    assert np.min(arr) == pytest.approx(0.0)
    assert np.sum(arr_counts) == pytest.approx(28547112)
    assert np.max(arr_counts) == pytest.approx(2258)
    assert np.min(arr_counts) == pytest.approx(0.0)

    record_all_counts = True
    grid = flex.double(flex.grid(grid_size, grid_size, grid_size), 0)
    counts = flex.int(flex.grid(grid_size, grid_size, grid_size), 0)
    process_tof_experiment_list(
        experiments, max_resolution, grid, counts, nproc, record_all_counts
    )
    normalize_voxels(grid, counts)
    arr = flumpy.to_numpy(grid)
    arr_counts = flumpy.to_numpy(counts)
    assert np.sum(arr) == pytest.approx(88440.4449734628)
    assert np.max(arr) == pytest.approx(568.0)
    assert np.min(arr) == pytest.approx(0.0)
    assert np.sum(arr_counts) == pytest.approx(33258685)
    assert np.max(arr_counts) == pytest.approx(2258)
    assert np.min(arr_counts) == pytest.approx(0.0)
