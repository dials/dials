from __future__ import annotations

from dxtbx.serialize import load
from scitbx import matrix

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
    expected = matrix.col((-0.627963, -0.778243, 0))

    assert axis.angle(expected) < 1e-6
