import procrunner


def test_dose_analysis_dials_data(dials_data, tmpdir):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # set cutoff to 0.0 to force one to be 'rejected'
    command = ["dials.dose_analysis", refls, expts, "min_completeness=0.4"]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.dose_analysis.html").check()


def test_dose_analysis_mtz(dials_data, tmpdir):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # First export the data
    command = ["dials.export", refls, expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.mtz").check()

    # set cutoff to 0.0 to force one to be 'rejected'
    command = [
        "dials.dose_analysis",
        "mtzfile=%s" % tmpdir.join("scaled.mtz").strpath,
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.dose_analysis.html").check()


def test_dose_analysis_input_handling(dials_data, tmpdir):
    """Test that errors are handled if more than one refl file, no refl/expt
    file or unscaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    command = ["dials.dose_analysis", refls, expts, refls]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr

    command = ["dials.dose_analysis", expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr

    command = ["dials.dose_analysis", refls]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr

    data_dir = dials_data("l_cysteine_dials_output")
    refls = data_dir / "20_integrated.pickle"
    expts = data_dir / "20_integrated_experiments.json"

    command = ["dials.dose_analysis", refls, expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode and result.stderr
