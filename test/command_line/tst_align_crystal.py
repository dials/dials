from __future__ import absolute_import, division

def run():
  import os
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.test_utils import show_diff
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(dials_regression, "experiment_test_data")
  cmd = "dials.align_crystal %s/kappa_experiments.json" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  out = "\n".join(result.stdout_lines[6:])
  assert not show_diff(out, """\
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
""")


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
