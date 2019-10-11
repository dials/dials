from __future__ import absolute_import, division, print_function

import os

import procrunner
import pytest
from cctbx import sgtbx

from dials.array_family import flex


def test_symmetry(dials_regression, tmpdir):
    """Simple test to check that dials.symmetry completes"""

    result = procrunner.run(
        [
            "dials.symmetry",
            os.path.join(dials_regression, "xia2-28", "20_integrated_experiments.json"),
            os.path.join(dials_regression, "xia2-28", "20_integrated.pickle"),
            os.path.join(dials_regression, "xia2-28", "25_integrated_experiments.json"),
            os.path.join(dials_regression, "xia2-28", "25_integrated.pickle"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check()
    assert tmpdir.join("symmetrized.expt").check()


from dials.algorithms.symmetry.cosym._generate_test_data import (
    generate_experiments_reflections,
)


def test_symmetry_R3m_I2m(tmpdir):
    os.chdir(tmpdir.strpath)
    space_group = "C 2"
    unit_cell = (53.173, 61.245, 69.292, 90.0, 93.04675, 90.0)
    space_group = sgtbx.space_group_info(space_group).group()
    experiments, reflections, _ = generate_experiments_reflections(
        space_group=space_group,
        unit_cell=unit_cell,
        sample_size=1,
        map_to_minimum=False,
    )
    experiments.as_json("tmp.expt")
    expt_file = tmpdir.join("tmp.expt").strpath
    joint_table = flex.reflection_table()
    for r in reflections:
        joint_table.extend(r)
    joint_table.as_pickle("tmp.refl")
    refl_file = tmpdir.join("tmp.refl").strpath

    # okay so generated in C 1 m 1, but lattice looks like R-3m, s*o reindex
    # to this first

    command = ["dials.symmetry", expt_file, refl_file]
    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)
    from dxtbx.serialize import load

    expts = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    for v, expected in zip(expts[0].crystal.get_unit_cell().parameters(), unit_cell):
        assert v == pytest.approx(expected)
