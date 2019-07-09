from __future__ import absolute_import, division, print_function

import pytest
import os
from libtbx import easy_run
from dxtbx.serialize import load


def test_modify_geometry(dials_regression, run_in_tmpdir):
    orig_expt_json = os.path.join(
        dials_regression, "experiment_test_data/kappa_experiments.json"
    )

    new_expt_json = os.path.join(os.getcwd(), "modified.expt")

    cmd = "dials.modify_geometry %s angles=10,20,30" % orig_expt_json
    result = easy_run.fully_buffered(cmd).raise_if_errors()

    assert os.path.exists(orig_expt_json), orig_expt_json
    orig_expt = load.experiment_list(orig_expt_json, check_format=False)
    assert os.path.exists(new_expt_json), new_expt_json
    new_expt = load.experiment_list(new_expt_json, check_format=False)

    orig_gonio = orig_expt.goniometers()[0]
    new_gonio = new_expt.goniometers()[0]
    assert orig_gonio.get_angles() == pytest.approx([0, 180, 0])
    assert new_gonio.get_angles() == pytest.approx([10, 20, 30])
