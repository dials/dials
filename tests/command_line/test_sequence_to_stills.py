from __future__ import annotations

import shutil
import subprocess

from dxtbx.model.experiment_list import ExperimentListFactory


def test_sequence_to_stills(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed", pathlib=True)
    input_experiments = data_dir / "integrated.expt"
    input_reflections = data_dir / "refined.refl"
    result = subprocess.run(
        [
            shutil.which("dials.sequence_to_stills"),
            input_experiments,
            input_reflections,
            "domain_size_ang=500",
            "half_mosaicity_deg=0.1",
            "max_scan_points=10",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "stills.expt").is_file()
    assert (tmp_path / "stills.refl").is_file()

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "stills.expt", check_format=False
    )
    assert len(experiments) == 10
    assert len(experiments.identifiers()) == 10


def test_data_with_static_model(dials_data, tmp_path):
    # test for regression of https://github.com/dials/dials/issues/2516
    data_dir = dials_data("insulin_processed", pathlib=True)
    input_experiments = data_dir / "indexed.expt"
    input_reflections = data_dir / "indexed.refl"
    result = subprocess.run(
        [
            shutil.which("dials.sequence_to_stills"),
            input_experiments,
            input_reflections,
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "stills.expt").is_file()
    assert (tmp_path / "stills.refl").is_file()

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "stills.expt", check_format=False
    )
    assert len(experiments) == 45


def test_sliced_sequence(dials_data, tmp_path):
    # test for regression of https://github.com/dials/dials/issues/2519
    data_dir = dials_data("insulin_processed", pathlib=True)
    input_experiments = data_dir / "indexed.expt"
    input_reflections = data_dir / "indexed.refl"

    result = subprocess.run(
        [
            shutil.which("dials.slice_sequence"),
            input_experiments,
            input_reflections,
            "image_range=5,45",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [
            shutil.which("dials.sequence_to_stills"),
            "indexed_5_45.expt",
            "indexed_5_45.refl",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "stills.expt").is_file()
    assert (tmp_path / "stills.refl").is_file()

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "stills.expt", check_format=False
    )
    assert len(experiments) == 41
