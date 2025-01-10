from __future__ import annotations

from libtbx import easy_run


def test_compare_orientation_matrices(dials_data, run_in_tmp_path):
    data_dir = dials_data("refinement_test_data", pathlib=True)
    cmd = (
        "dials.compare_orientation_matrices %s/i04-weak.json %s/i04-weak-regression.json"
        % (data_dir, data_dir)
    )
    result = easy_run.fully_buffered(cmd).raise_if_errors()
    out = "\n".join(result.stdout_lines[7:])
    out = out.replace("-0", "0")
    assert (
        out
        == """\
Change of basis op: a,b,c
Rotation matrix to transform crystal 1 to crystal 2:
{{1.000, 0.000, 0.000},
 {0.000, 1.000, 0.000},
 {0.000, 0.000, 1.000}}
Rotation of 0.002 degrees about axis (0.917, 0.082, 0.390)
"""
    )
