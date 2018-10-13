from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run

def test_compare_orientation_matrices(dials_regression, run_in_tmpdir):
  path = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  cmd = "dials.compare_orientation_matrices %s/experiments.json %s/regression_experiments.json" %(path, path)
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  out = "\n".join(result.stdout_lines[7:])
  out = out.replace("-0", "0")
  assert out == """\
Change of basis op: a,b,c
Rotation matrix to transform crystal 1 to crystal 2:
{{1.000, 0.000, 0.000},
 {0.000, 1.000, 0.000},
 {0.000, 0.000, 1.000}}
Rotation of 0.002 degrees about axis (0.916, 0.081, 0.393)
"""
