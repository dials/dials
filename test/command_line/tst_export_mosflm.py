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

  path = os.path.join(dials_regression, "experiment_test_data")

  cmd = "dials.export_mosflm '%s/experiment_1.json'" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("mosflm/index.mat")
  with open("mosflm/index.mat", "rb") as f:
    lines = f.read()
    assert not show_diff(lines, """\
 -0.01210200 -0.01954526  0.00309519
 -0.00416605 -0.00080573 -0.02427340
  0.01931593 -0.01241956 -0.00329641
       0.000       0.000       0.000
 -0.57227172 -0.36749719  0.25435065
 -0.19700176  0.09729719 -0.77371631
  0.91339955 -1.06160155 -1.71101839
     42.2717     42.2720     39.6704     90.0001     89.9993     89.9998
       0.000       0.000       0.000
""")
  assert os.path.exists("mosflm/mosflm.in")
  with open("mosflm/mosflm.in", "rb") as f:
    lines = f.read()
    assert lines.startswith("DIRECTORY")
    assert not show_diff("\n".join(lines.split("\n")[1:]), """\
TEMPLATE centroid_####.cbf
SYMMETRY 89
BEAM 220.002 212.478
DISTANCE 190.1800
MATRIX index.mat
""")


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
