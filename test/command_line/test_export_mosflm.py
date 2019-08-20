from __future__ import absolute_import, division, print_function

import json
import os

import procrunner


def test_export_mosflm(dials_regression, tmpdir):
    dials_regression_escaped = json.dumps(dials_regression).strip('"')
    with open(
        os.path.join(dials_regression, "experiment_test_data/experiment_1.json"), "r"
    ) as fi:
        with (tmpdir / "experiments.json").open("w") as fo:
            fo.write(fi.read().replace("$DIALS_REGRESSION", dials_regression_escaped))

    tmpdir.chdir()

    result = procrunner.run(["dials.export", "format=mosflm", "experiments.json"])
    assert not result.returncode and not result.stderr

    assert os.path.exists("mosflm/index.mat")
    with open("mosflm/index.mat", "rb") as f:
        lines = f.read()
    assert (
        lines
        == """
 -0.01210200 -0.01954526  0.00309519
 -0.00416605 -0.00080573 -0.02427340
  0.01931593 -0.01241956 -0.00329641
       0.000       0.000       0.000
 -0.52228050 -0.84350975  0.12535704
 -0.17980379 -0.03477015 -0.98308781
  0.83360283 -0.53598726 -0.13350648
     42.2717     42.2720     39.6704     90.0001     89.9993     89.9998
       0.000       0.000       0.000
""".strip(
            "\n"
        )
    )
    assert os.path.exists("mosflm/mosflm.in")
    with open("mosflm/mosflm.in", "rb") as f:
        lines = f.read()
    assert (
        lines
        == """
DIRECTORY %s%scentroid_test_data
TEMPLATE centroid_####.cbf
SYMMETRY 89
BEAM 220.002 212.478
DISTANCE 190.1800
MATRIX index.mat
""".strip(
            "\n"
        )
        % (dials_regression, os.path.sep)
    )
