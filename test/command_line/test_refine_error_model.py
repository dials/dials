"""
Test standalone error model refinement on scaled data
"""
import procrunner


def test_refine_error_model(dials_data, tmpdir):
    """Test the program on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl")
    expts = location.join("scaled_20_25.expt")

    # Use a range of options, some default.
    command = [
        "dials.refine_error_model",
        refls,
        expts,
        "json=error_model.json",
        "html=error_model.html",
        "intensity_choice=combine",
        "combine.Imid=250",
        "basic.minimisation=individual",
    ]

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("error_model.html").check(file=1)
    assert tmpdir.join("error_model.json").check(file=1)
