from __future__ import annotations

from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.serialize import load
from scitbx import matrix
from scitbx.array_family import flex

import dials.command_line.find_rotation_axis as dials_find_rotation_axis


def test_find_rotation_axis(dials_data, run_in_tmp_path):
    myd88 = dials_data("MyD88_processed", pathlib=True)
    dials_find_rotation_axis.run(
        args=[str(myd88 / "imported.expt"), str(myd88 / "strong.refl")]
    )
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()
    assert run_in_tmp_path.joinpath("find_rotation_axis-histogram.png").is_file()
    assert run_in_tmp_path.joinpath("find_rotation_axis-projection.png").is_file()
    experiments = load.experiment_list(
        run_in_tmp_path / "optimised.expt", check_format=False
    )

    axis = matrix.col(experiments[0].goniometer.get_rotation_axis())
    expected = matrix.col((-0.626604, -0.779338, 0))

    assert axis.angle(expected) < 1e-6


def test_find_rotation_axis_multi_axis_goniometer(dials_data, run_in_tmp_path):
    myd88 = dials_data("MyD88_processed", pathlib=True)
    experiments = load.experiment_list(str(myd88 / "imported.expt"), check_format=False)

    # Replace single-axis goniometer with multi-axis goniometer
    axes = flex.vec3_double([(0, -1, 0), (0, -0.642788, 0.766044), (0, -1, 0)])
    angles = flex.double((0, 0, -62))
    names = flex.std_string(("PHI", "KAPPA=CHI", "OMEGA"))
    scan_axis = 2
    experiments[0].goniometer = GoniometerFactory.multi_axis(
        axes, angles, names, scan_axis
    )
    experiments.as_file("modified.expt")

    # Now run find_rotation_axis, and expect the same result.
    dials_find_rotation_axis.run(args=["modified.expt", str(myd88 / "strong.refl")])
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()
    assert run_in_tmp_path.joinpath("find_rotation_axis-histogram.png").is_file()
    assert run_in_tmp_path.joinpath("find_rotation_axis-projection.png").is_file()
    experiments = load.experiment_list(
        run_in_tmp_path / "optimised.expt", check_format=False
    )

    axis = matrix.col(experiments[0].goniometer.get_rotation_axis())
    expected = matrix.col((-0.626604, -0.779338, 0))

    assert axis.angle(expected) < 1e-6
