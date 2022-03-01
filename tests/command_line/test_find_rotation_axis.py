from __future__ import annotations

import os

from dxtbx.serialize import load
from scitbx import matrix

import dials.command_line.find_rotation_axis as dials_find_rotation_axis


def test_find_rotation_axis(dials_data, run_in_tmpdir):
    myd88 = dials_data("MyD88_processed")
    args = [myd88.join("imported.expt").strpath, myd88.join("strong.refl").strpath]
    dials_find_rotation_axis.run(args=args)
    assert os.path.isfile("optimised.expt")
    assert os.path.isfile("find_rotation_axis-histogram.png")
    assert os.path.isfile("find_rotation_axis-projection.png")
    experiments = load.experiment_list("optimised.expt", check_format=False)

    axis = matrix.col(experiments[0].goniometer.get_rotation_axis())
    expected = matrix.col((-0.627963, -0.778243, 0))

    assert axis.angle(expected) < 1e-6
