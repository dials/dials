"""
Test command line program dials.two_theta_refine by running a job with saved
data and comparing with expected output.
"""


from __future__ import annotations

from pathlib import Path

import procrunner
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory


def test(dials_data, tmpdir):
    """Test two theta refine on integrated data."""
    # use multiple scan small molecule data for this test
    data_dir = dials_data("l_cysteine_dials_output", pathlib=True)
    prefix = (20, 25, 30, 35)
    exp_path = [data_dir / ("%d_integrated_experiments.json" % p) for p in prefix]
    pkl_path = [data_dir / ("%d_integrated.pickle" % p) for p in prefix]

    for pth in exp_path + pkl_path:
        assert pth.is_file(), f"{str(pth)} missing"

    cmd = (
        [
            "dials.two_theta_refine",
            "cif=refined_cell.cif",
            "output.correlation_plot.filename=corrplot.png",
        ]
        + [str(p) for p in exp_path]
        + [str(p) for p in pkl_path]
    )

    print(cmd)

    # work in a temporary directory
    result = procrunner.run(cmd, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert Path(tmpdir / "refined_cell.expt").is_file()
    ref_exp = ExperimentListFactory.from_json_file(
        str(tmpdir / "refined_cell.expt"), check_format=False
    )

    xls = ref_exp.crystals()
    assert len(xls) == 4
    for xl in xls:
        assert xl.get_unit_cell() != xl.get_recalculated_unit_cell()
        # test refined crystal model against expected values
        assert xl.get_recalculated_unit_cell().parameters() == pytest.approx(
            (5.428022880, 8.144145476, 12.039666971, 90.0, 90.0, 90.0), 1e-4
        )
        assert xl.get_recalculated_cell_parameter_sd() == pytest.approx(
            (9.58081e-5, 0.000149909, 0.000215765, 0, 0, 0), 1e-4
        )
        assert xl.get_recalculated_cell_volume_sd() == pytest.approx(0.0116254298, 1e-4)


def test_two_theta_refine_scaled_data(dials_data, tmpdir):
    """Test two theta refine on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = str(location / "scaled_20_25.refl")
    expts = str(location / "scaled_20_25.expt")

    command = [
        "dials.two_theta_refine",
        refls,
        expts,
        "output.experiments=refined_cell.expt",
        "partiality_threshold=0.99",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert Path(tmpdir / "refined_cell.expt").is_file()

    ref_exp = ExperimentListFactory.from_json_file(
        str(tmpdir / "refined_cell.expt"), check_format=False
    )

    assert len(ref_exp.crystals()) == 2
    for xl in ref_exp.crystals():
        # test refined crystal model against expected values
        assert xl.get_recalculated_unit_cell().parameters() == pytest.approx(
            (5.426921, 8.146654, 12.037366, 90.0, 90.0, 90.0), 1e-4
        )
        assert xl.get_recalculated_cell_parameter_sd() == pytest.approx(
            (2.0123e-04, 2.8039e-04, 4.5284e-04, 0, 0, 0), 1e-4
        )
        assert xl.get_recalculated_cell_volume_sd() == pytest.approx(0.0237364, 1e-4)
