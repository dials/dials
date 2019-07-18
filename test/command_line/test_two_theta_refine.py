"""
Test command line program dials.two_theta_refine by running a job with saved
data and comparing with expected output.
"""

from __future__ import absolute_import, division, print_function

import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentListFactory


def test(dials_regression, run_in_tmpdir):
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
        "dials.two_theta_refine "
        + " ".join(exp_path)
        + " "
        + " ".join(pkl_path)
        + " cif=refined_cell.cif"
        + " "
        "output.correlation_plot.filename=corrplot.png"
    )

    print(cmd)

    # work in a temporary directory
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    ref_exp = ExperimentListFactory.from_json_file(
        "refined_cell.expt", check_format=False
    )

    xls = ref_exp.crystals()
    assert len(xls) == 1  # crystal models should have been combined
    xl = xls[0]

    # test refined crystal model against expected values
    assert approx_equal(
        xl.get_unit_cell().parameters(),
        (5.428022880, 8.144145476, 12.039666971, 90.0, 90.0, 90.0),
    )
    assert approx_equal(
        xl.get_cell_parameter_sd(), (9.58081e-5, 0.000149909, 0.000215765, 0, 0, 0)
    )
    assert approx_equal(xl.get_cell_volume_sd(), 0.0116254298)
