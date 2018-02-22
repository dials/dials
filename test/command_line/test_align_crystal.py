from __future__ import absolute_import, division

import pytest
import os
from libtbx import easy_run

def test_align_crystal(dials_regression, tmpdir):
  tmpdir.chdir()

  path = os.path.join(dials_regression, "experiment_test_data")
  cmd = "dials.align_crystal %s/kappa_experiments.json" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  out = "\n".join(result.stdout_lines[6:])
  assert out == """\
Angles between reciprocal cell axes and principal experimental axes:
--------------------------------------------
Experimental axis | a*     | b*     | c*
--------------------------------------------
GON_PHI           | 88.712 | 88.317 | 1.760
Beam              | 13.904 | 46.197 | 88.481
GON_OMEGA         | 88.712 | 88.317 | 1.760
--------------------------------------------

Angles between unit cell axes and principal experimental axes:
--------------------------------------------
Experimental axis | a      | b      | c
--------------------------------------------
GON_PHI           | 89.484 | 88.801 | 1.760
Beam              | 43.844 | 76.182 | 88.481
GON_OMEGA         | 89.484 | 88.801 | 1.760
--------------------------------------------

Independent solutions:
----------------------------------------------------
Primary axis | Secondary axis | GON_KAPPA | GON_PHI
----------------------------------------------------
c* (6-fold)  | a* (2-fold)    |   4.324   | -77.874
c* (6-fold)  | a* (2-fold)    |  -4.324   |  106.075
----------------------------------------------------\
"""
