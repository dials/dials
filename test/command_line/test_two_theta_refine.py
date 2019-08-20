"""
Test command line program dials.two_theta_refine by running a job with saved
data and comparing with expected output.
"""

from __future__ import absolute_import, division, print_function

import os
import procrunner
import pytest
from dxtbx.model.experiment_list import ExperimentListFactory


def test(dials_regression, tmpdir):
    """Test two theta refine on integrated data."""
    # use multiple scan small molecule data for this test
    data_dir = os.path.join(dials_regression, "xia2-28")
    prefix = ["20", "25", "30", "35"]
    exp_path = [e + "_integrated_experiments.json" for e in prefix]
    exp_path = [os.path.join(data_dir, e) for e in exp_path]
    pkl_path = [e + "_integrated.pickle" for e in prefix]
    pkl_path = [os.path.join(data_dir, e) for e in pkl_path]

    for pth in exp_path + pkl_path:
        assert os.path.exists(pth), "%s missing" % pth

    cmd = (
        [
            "dials.two_theta_refine",
            "cif=refined_cell.cif",
            "output.correlation_plot.filename=corrplot.png",
        ]
        + exp_path
        + pkl_path
    )

    print(cmd)

    # work in a temporary directory
    result = procrunner.run(cmd, working_directory=tmpdir)
    assert result.returncode == 0
    assert result.stderr == ""
    assert tmpdir.join("refined_cell.expt").check()
    ref_exp = ExperimentListFactory.from_json_file(
        tmpdir.join("refined_cell.expt").strpath, check_format=False
    )

    xls = ref_exp.crystals()
    assert len(xls) == 1  # crystal models should have been combined
    xl = xls[0]

    # test refined crystal model against expected values
    assert xl.get_unit_cell().parameters() == pytest.approx(
        (5.428022880, 8.144145476, 12.039666971, 90.0, 90.0, 90.0), 1e-4
    )
    assert xl.get_cell_parameter_sd() == pytest.approx(
        (9.58081e-5, 0.000149909, 0.000215765, 0, 0, 0), 1e-4
    )
    assert xl.get_cell_volume_sd() == pytest.approx(0.0116254298, 1e-4)


def test_two_theta_refine_scaled_data(dials_data, tmpdir):
    """Test two theta refine on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")

    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    command = [
        "dials.two_theta_refine",
        refls,
        expts,
        "output.experiments=refined_cell.expt",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode == 0
    assert result.stderr == ""
    assert tmpdir.join("refined_cell.expt").check()

    ref_exp = ExperimentListFactory.from_json_file(
        tmpdir.join("refined_cell.expt").strpath, check_format=False
    )

    assert len(ref_exp.crystals()) == 1  # crystal models should have been combined
    xl = ref_exp.crystals()[0]

    # test refined crystal model against expected values
    assert xl.get_unit_cell().parameters() == pytest.approx(
        (5.426921, 8.146654, 12.037366, 90.0, 90.0, 90.0), 1e-4
    )
    assert xl.get_cell_parameter_sd() == pytest.approx(
        (2.0123e-04, 2.8039e-04, 4.5284e-04, 0, 0, 0), 1e-4
    )
    assert xl.get_cell_volume_sd() == pytest.approx(0.0237364, 1e-4)
