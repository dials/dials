from __future__ import division

def run():

  import os
  import shutil
  import fileinput
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.test_utils import approx_equal
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)

  orig_expt_json = os.path.join(
    dials_regression, "experiment_test_data/kappa_experiments.json")

  new_expt_json = os.path.join(os.getcwd(), 'modified_experiments.json')

  cmd = "dials.modify_geometry %s angles=10,20,30" %orig_expt_json
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  from dxtbx.serialize import load
  assert os.path.exists(orig_expt_json), orig_expt_json
  orig_expt = load.experiment_list(orig_expt_json, check_format=False)
  assert os.path.exists(new_expt_json), new_expt_json
  new_expt = load.experiment_list(new_expt_json, check_format=False)

  orig_gonio = orig_expt.goniometers()[0]
  new_gonio = new_expt.goniometers()[0]
  assert approx_equal(orig_gonio.get_angles(), [0,180,0])
  assert approx_equal(new_gonio.get_angles(), [10,20,30])

  return


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
