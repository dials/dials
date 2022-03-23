from __future__ import annotations

from os import path

import procrunner
import pytest

from dxtbx.serialize import load

from dials.command_line.modify_geometry import update


def test_run(dials_regression, tmp_path):
    orig_expt_json = path.join(
        dials_regression, "experiment_test_data", "kappa_experiments.json"
    )
    assert path.exists(orig_expt_json)

    orig_expt = load.experiment_list(orig_expt_json, check_format=False)
    orig_gonio = orig_expt.goniometers()[0]
    assert orig_gonio.get_angles() == pytest.approx([0, 180, 0])

    result = procrunner.run(
        ["dials.modify_geometry", orig_expt_json, "angles=10,20,30"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    new_expt_json = tmp_path / "modified.expt"
    assert new_expt_json.is_file()

    new_expt = load.experiment_list(new_expt_json, check_format=False)

    new_gonio = new_expt.goniometers()[0]
    assert new_gonio.get_angles() == pytest.approx([10, 20, 30])


def test_load(dials_data, tmp_path):
    from dials.command_line.modify_geometry import phil_scope as modify_geom_phil

    orig_expt = dials_data("aluminium_standard", pathlib=True) / "imported.expt"
    assert orig_expt.is_file()

    orig_expt = load.experiment_list(orig_expt, check_format=False)
    orig_beam = orig_expt.beams()[0]
    assert orig_beam.get_wavelength() == pytest.approx(0.02508235604)

    new_expt = tmp_path / "modified.expt"

    ## this is the idiosyncratic way to update phils
    ## and it should work but in fact does nothing
    working_params = modify_geom_phil.extract()
    working_params.geometry.beam.wavelength = 0.05
    working_params.output.experiments = str(new_expt)

    update(orig_expt, working_params)
    assert new_expt.is_file()

    new_expt = load.experiment_list(new_expt, check_format=False)

    new_beam = new_expt.beams()[0]
    assert new_beam.get_wavelength() == pytest.approx(0.05)
