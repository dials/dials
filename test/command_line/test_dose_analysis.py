"""Tests for dials.dose_analysis"""
import os
import procrunner
import pytest
from dials.array_family import flex
from dials.command_line.dose_analysis import PychefRunner, phil_scope, run
from dxtbx.serialize import load


def test_dose_analysis_dials_data(dials_data, run_in_tmpdir):
    """Test dials.dose_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    args = [
        refls,
        expts,
        "min_completeness=0.4",
        "-v",
        "json=dials.dose_analysis.json",
    ]
    run(args)
    assert os.path.isfile("dials.dose_analysis.html")
    assert os.path.isfile("dials.dose_analysis.json")


def test_setup_from_dials_data(dials_data, run_in_tmpdir):
    """Test dials.dose_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath
    table = flex.reflection_table.from_file(refls)
    experiments = load.experiment_list(expts, check_format=False)
    params = phil_scope.extract()
    params.dose.experiments.shared_crystal = True
    params.dose.experiments.dose_per_image = [1.0, 2.0]
    # First experiment is images 1>1800, second 1>1700 (i.e. dose spans 1>5200)
    runner = PychefRunner.from_dials_datafiles(params, experiments, table)
    assert max(runner.dose) == 5198  # last reflection measured not quite at end
    assert min(runner.dose) == 2

    # Now try again in 'standard' mode i.e. not shared crystal, and set a
    # starting dose
    params.dose.experiments.shared_crystal = False
    params.dose.experiments.dose_per_image = [1.0]
    params.dose.experiments.starting_doses = [10, 10]
    runner = PychefRunner.from_dials_datafiles(params, experiments, table)
    assert max(runner.dose) == 1800 + 10
    assert min(runner.dose) == 2 + 10


def test_dose_analysis_mtz(dials_data, run_in_tmpdir):
    """Test dials.dose_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # First export the data
    command = ["dials.export", refls, expts]
    result = procrunner.run(command)
    assert not result.returncode and not result.stderr
    assert os.path.isfile("scaled.mtz")

    args = [
        run_in_tmpdir.join("scaled.mtz").strpath,
        "anomalous=True",
        "json=dials.dose_analysis.json",
    ]
    run(args)
    assert os.path.isfile("dials.dose_analysis.html")
    assert os.path.isfile("dials.dose_analysis.json")


def test_dose_analysis_input_handling(dials_data, run_in_tmpdir):
    """Test that errors are handled if more than one refl file, no refl/expt
    file or unscaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

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

    # Unscaled data
    location = dials_data("l_cysteine_dials_output")
    refls = location.join("20_integrated.pickle").strpath
    expts = location.join("20_integrated_experiments.json").strpath

    args = [refls, expts]
    with pytest.raises(SystemExit):
        run(args)
