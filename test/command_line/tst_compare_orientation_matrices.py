from __future__ import division

def run():
  import os
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.test_utils import show_diff
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  cmd = "dials.compare_orientation_matrices '%s/experiments.json' '%s/regression_experiments.json'" %(path, path)
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  expected_out = [
    'Rotation matrix to transform crystal 1 to crystal 2',
    '{{0.999999998813, -3.68587676544*^-05, 3.18565817573*^-05},',
    ' {3.6859180576*^-05, 0.999999999237, -1.29613953456*^-05},',
    ' {-3.18561039918*^-05, 1.2962569538*^-05, 0.999999999409}}',
    'Euler angles (xyz): 0.00, 0.00, 0.00', ''
  ]
  for i, line in enumerate(result.stdout_lines[-6:]):
    assert not show_diff(line, expected_out[i])


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
