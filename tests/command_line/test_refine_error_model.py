from __future__ import annotations

import procrunner


def test_standalone_error_model_refinement_on_scaled_data(dials_data, tmp_path):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

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

    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("error_model.html").is_file()
    assert tmp_path.joinpath("error_model.json").is_file()
