from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from cctbx import crystal, sgtbx, uctbx
import pytest


def test_symmetry_analysis():
    coords = flex.double(
        [
            [0.835, 0.158],
            [0.772, 0.104],
            [0.108, 0.907],
            [0.058, 0.76],
            [0.926, 0.189],
            [0.221, 0.888],
            [0.957, 0.137],
            [0.958, 0.143],
            [-0.015, 0.726],
            [-0.066, 0.29],
            [0.135, 0.848],
            [0.085, 0.788],
            [0.897, 0.126],
            [0.749, 0.073],
            [0.166, 0.943],
            [0.871, 0.248],
            [0.116, 0.968],
            [0.116, 0.973],
            [0.706, 0.007],
            [0.288, -0.055],
            [0.137, 0.848],
            [0.089, 0.78],
            [0.893, 0.122],
            [0.749, 0.077],
            [0.165, 0.941],
            [0.877, 0.242],
            [0.114, 0.968],
            [0.12, 0.971],
            [0.716, 0.002],
            [0.292, -0.062],
            [0.841, 0.162],
            [0.774, 0.104],
            [0.1, 0.909],
            [0.054, 0.761],
            [0.927, 0.184],
            [0.227, 0.88],
            [0.957, 0.137],
            [0.961, 0.143],
            [-0.007, 0.716],
            [-0.061, 0.287],
            [0.13, 0.848],
            [0.084, 0.783],
            [0.898, 0.124],
            [0.749, 0.075],
            [0.169, 0.94],
            [0.871, 0.247],
            [0.114, 0.969],
            [0.12, 0.969],
            [0.717, 0.0],
            [0.296, -0.066],
            [0.84, 0.154],
            [0.776, 0.103],
            [0.104, 0.908],
            [0.057, 0.755],
            [0.925, 0.19],
            [0.227, 0.883],
            [0.958, 0.136],
            [0.962, 0.143],
            [-0.017, 0.724],
            [-0.067, 0.295],
        ]
    )

    sym_ops = [
        sgtbx.rt_mx(s)
        for s in ("-z,-y,-x", "y,z,x", "x,y,z", "-x,-z,-y", "z,x,y", "-y,-x,-z")
    ]

    crystal_symmetry = crystal.symmetry(
        unit_cell=uctbx.unit_cell((98.33, 98.33, 135.99, 90, 90, 120)),
        space_group_info=sgtbx.space_group_info("R3:H"),
    ).minimum_cell()

    from cctbx.sgtbx.lattice_symmetry import metric_subgroups

    subgroups = metric_subgroups(
        crystal_symmetry, max_delta=5, bravais_types_only=False
    )

    cb_op_inp_min = sgtbx.change_of_basis_op()

    from dials.algorithms.symmetry.cosym import SymmetryAnalysis

    analysis = SymmetryAnalysis(coords, sym_ops, subgroups, cb_op_inp_min)

    assert analysis.best_solution.likelihood > 0.99
    assert analysis.best_solution.confidence > 0.98
    assert (
        analysis.best_solution.subgroup["best_subsym"].space_group().type().number()
        == 148
    )  # R -3 :H
    assert (
        str(analysis)
        == """\
Scoring individual symmetry elements
----------------------------------------------
likelihood  Z-CC   CC         Operator
----------------------------------------------
0.087       1.96   0.20        2 |(0, -1, 1)
0.087       1.96   0.20        2 |(-1, 0, 1)
0.949       10.00  1.00  ***   3^-1 |(1, 1, 1)
0.087       1.96   0.20        2 |(-1, 1, 0)
0.949       10.00  1.00  ***   3 |(1, 1, 1)
----------------------------------------------
Scoring all possible sub-groups
--------------------------------------------------------------------------------
Patterson group       Likelihood  NetZcc  Zcc+    Zcc-   delta  Reindex operator
--------------------------------------------------------------------------------
R -3 :H          ***  0.995        8.04    10.00   1.96  0.0    b-c,-a+c,a+b+c
P -1                  0.003       -6.50    0.00    6.50  0.0    a,b,c
R -3 m :H             0.001        6.50    6.50    0.00  0.0    b-c,-a+c,a+b+c
C 1 2/m 1             0.000       -5.24    1.96    7.21  0.0    -a-b,a-b,c
C 1 2/m 1             0.000       -5.24    1.96    7.21  0.0    -b-c,b-c,a
C 1 2/m 1             0.000       -5.24    1.96    7.21  0.0    -a-c,-a+c,b
--------------------------------------------------------------------------------
Best solution: R -3 :H
Unit cell: (98.33, 98.33, 135.99, 90, 90, 120)
Reindex operator: b-c,-a+c,a+b+c
Laue group probability: 0.995
Laue group confidence: 0.994"""
    )

    d = analysis.as_dict()
    assert d["sym_op_scores"][0] == {
        "cc": pytest.approx(0.19620531091685714),
        "operator": "-x,-z,-y",
        "likelihood": pytest.approx(0.08665625555575088),
        "stars": "",
        "z_cc": pytest.approx(1.9620531091685713),
    }
    assert d["subgroup_scores"][0] == {
        "confidence": pytest.approx(0.9940687431995551),
        "z_cc_for": pytest.approx(9.999725360190128),
        "stars": "***",
        "patterson_group": "-R 3",
        "max_angular_difference": 0.0,
        "likelihood": pytest.approx(0.995493024305035),
        "cb_op": "-1/3*x+2/3*y-1/3*z,-2/3*x+1/3*y+1/3*z,1/3*x+1/3*y+1/3*z",
        "z_cc_against": pytest.approx(1.9620621986200772),
        "unit_cell": pytest.approx(
            (
                98.32999999999998,
                98.32999999999998,
                135.99,
                90.0,
                90.0,
                119.99999999999999,
            )
        ),
        "z_cc_net": pytest.approx(8.037663161570052),
    }
