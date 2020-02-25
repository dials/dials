"""Tests for dials.dose_analysis"""
import procrunner
from dials.command_line.dose_analysis import PychefRunner, phil_scope
from dxtbx.serialize import load
from dials.array_family import flex


def test_dose_analysis_dials_data(dials_data, tmpdir):
    """Test dials.dose_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    command = [
        "dials.dose_analysis",
        refls,
        expts,
        "min_completeness=0.4",
        "-v",
        "json=dials.dose_analysis.json",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.dose_analysis.html").check()
    assert tmpdir.join("dials.dose_analysis.json").check()


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


def test_dose_analysis_mtz(dials_data, tmpdir):
    """Test dials.dose_analysis on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # First export the data
    command = ["dials.export", refls, expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.mtz").check()

    command = [
        "dials.dose_analysis",
        tmpdir.join("scaled.mtz").strpath,
        "anomalous=True",
        "json=dials.dose_analysis.json",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.dose_analysis.html").check()
    assert tmpdir.join("dials.dose_analysis.json").check()


def test_dose_analysis_input_handling(dials_data, tmpdir):
    """Test that errors are handled if more than one refl file, no refl/expt
    file or unscaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # Too many refl files
    command = ["dials.dose_analysis", refls, expts, refls]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr

    # No refl file
    command = ["dials.dose_analysis", expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr

    # No expt file
    command = ["dials.dose_analysis", refls]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr

    # Unscaled data
    data_dir = dials_data("l_cysteine_dials_output")
    refls = data_dir / "20_integrated.pickle"
    expts = data_dir / "20_integrated_experiments.json"

    command = ["dials.dose_analysis", refls, expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr
