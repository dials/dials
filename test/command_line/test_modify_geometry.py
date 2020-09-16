from __future__ import absolute_import, division, print_function

import os

import procrunner
import pytest

from dxtbx.serialize import load


def test_modify_geometry(dials_regression, tmpdir):
    orig_expt_json = os.path.join(
        dials_regression, "experiment_test_data/kappa_experiments.json"
    )

    new_expt_json = tmpdir.join("modified.expt")

    result = procrunner.run(
        ["dials.modify_geometry", orig_expt_json, "angles=10,20,30"],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    assert os.path.exists(orig_expt_json), orig_expt_json
    orig_expt = load.experiment_list(orig_expt_json, check_format=False)
    assert new_expt_json.check()
    new_expt = load.experiment_list(new_expt_json.strpath, check_format=False)

    orig_gonio = orig_expt.goniometers()[0]
    new_gonio = new_expt.goniometers()[0]
    assert orig_gonio.get_angles() == pytest.approx([0, 180, 0])
    assert new_gonio.get_angles() == pytest.approx([10, 20, 30])
