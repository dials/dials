from __future__ import annotations

import os
import shutil
import subprocess

import pytest

from dxtbx.serialize import load

from dials.array_family import flex
from dials.command_line.damage_analysis import PychefRunner, phil_scope, run


def test_damage_analysis_on_scaled_data(dials_data, run_in_tmp_path):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = str(location / "scaled_20_25.refl")
    expts = str(location / "scaled_20_25.expt")

    args = [
        str(refls),
        str(expts),
        "min_completeness=0.4",
        "-v",
        "json=dials.damage_analysis.json",
    ]
    run(args)
    assert run_in_tmp_path.joinpath("dials.damage_analysis.html").is_file()
    assert run_in_tmp_path.joinpath("dials.damage_analysis.json").is_file()


def test_damage_analysis_damage_series(dials_data, run_in_tmp_path):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    args = [
        str(refls),
        str(expts),
        "min_completeness=0.4",
        "-v",
        "json=dials.damage_analysis.json",
        "dose_group_size=500",
        "output.damage_series=True",
        "output.accumulation_series=True",
    ]
    run(args)
    assert os.path.isfile("dials.damage_analysis.html")
    assert os.path.isfile("dials.damage_analysis.json")
    expected_series = [
        "0_500",
        "500_1000",
        "1000_1500",
        "1500_2000",
        "0_1000",
        "0_1500",
        "0_2000",
    ]
    for e in expected_series:
        assert os.path.isfile(f"damage_series_{e}.refl")
        assert os.path.isfile(f"damage_series_{e}.expt")


def test_setup_from_dials_data(dials_data):
    """Test dials.damage_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"
    table = flex.reflection_table.from_file(refls)
    experiments = load.experiment_list(expts, check_format=False)
    params = phil_scope.extract()
    params.dose.experiments.shared_crystal = True
    params.dose.experiments.dose_per_image = [1.0, 2.0]
    # First experiment is images 1>1800, second 1>1700 (i.e. dose spans 1>5200)
    runner = PychefRunner.from_dials_data_files(params, experiments, table)
    assert max(runner.dose) == 5198  # last reflection measured not quite at end
    assert min(runner.dose) == 2

    # Now try again in 'standard' mode i.e. not shared crystal, and set a
    # starting dose
    params.dose.experiments.shared_crystal = False
    params.dose.experiments.dose_per_image = [1.0]
    params.dose.experiments.starting_doses = [10, 10]
    runner = PychefRunner.from_dials_data_files(params, experiments, table)
    assert max(runner.dose) == 1800 + 10
    assert min(runner.dose) == 2 + 10


def test_damage_analysis_on_scaled_mtz(dials_data, run_in_tmp_path):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = str(location / "scaled_20_25.refl")
    expts = str(location / "scaled_20_25.expt")

    # First export the data
    command = [shutil.which("dials.export"), refls, expts]
    result = subprocess.run(command, cwd=run_in_tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert os.path.isfile("scaled.mtz")

    args = [
        str(run_in_tmp_path / "scaled.mtz"),
        "anomalous=True",
        "json=dials.damage_analysis.json",
    ]
    run(args)
    assert run_in_tmp_path.joinpath("dials.damage_analysis.html").is_file()
    assert run_in_tmp_path.joinpath("dials.damage_analysis.json").is_file()


def test_damage_analysis_mtz_damage_series(dials_data, run_in_tmp_path):
    """Test dials.damage_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    # First export the data
    command = [shutil.which("dials.export"), refls, expts]
    result = subprocess.run(command, cwd=run_in_tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert os.path.isfile("scaled.mtz")

    args = [
        str(run_in_tmp_path / "scaled.mtz"),
        "anomalous=True",
        "json=dials.damage_analysis.json",
        "dose_group_size=500",
        "output.damage_series=True",
        "output.accumulation_series=True",
    ]
    run(args)
    assert os.path.isfile("dials.damage_analysis.html")
    assert os.path.isfile("dials.damage_analysis.json")
    expected_series = [
        "0_500",
        "500_1000",
        "1000_1500",
        "1500_2000",
        "0_1000",
        "0_1500",
        "0_2000",
    ]
    for e in expected_series:
        assert os.path.isfile(f"damage_series_{e}.mtz")


def test_damage_analysis_input_handling(dials_data, run_in_tmp_path):
    """Test that errors are handled if more than one refl file, no refl/expt
    file or unscaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = str(location / "scaled_20_25.refl")
    expts = str(location / "scaled_20_25.expt")

    # Too many refl files
    args = [refls, expts, refls]
    with pytest.raises(SystemExit):
        run(args)

    # No refl file
    args = [expts]
    with pytest.raises(SystemExit):
        run(args)

    # No expt file
    args = [refls]
    with pytest.raises(SystemExit):
        run(args)


def test_damage_analysis_fails_on_unscaled_data(dials_data, run_in_tmp_path):
    location = dials_data("l_cysteine_dials_output", pathlib=True)
    refls = str(location / "20_integrated.pickle")
    expts = str(location / "20_integrated_experiments.json")

    args = [refls, expts]
    with pytest.raises(SystemExit):
        run(args)
