from __future__ import division

# this import required early to avoid seg fault on some systems
try:
  import scipy.linalg # import dependency
except ImportError, e:
  pass

import os
import libtbx.load_env
from libtbx.test_utils import approx_equal
from cctbx import uctbx

def run():
  have_dials_regression = libtbx.env.has_module("dials_regression")
  if not have_dials_regression:
    print "Skipped: dials_regression not available"
    return
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  from dials.test.algorithms.indexing.tst_index import run_one_indexing

  expected_unit_cell = uctbx.unit_cell(
    (11.624, 13.550, 30.103, 89.964, 93.721, 90.132))
  expected_rmsds = (0.039, 0.035, 0.002)

  datablock_old = os.path.join(
    dials_regression, "indexing_test_data/phi_scan/datablock_old.json")
  datablock_new = os.path.join(
    dials_regression, "indexing_test_data/phi_scan/datablock.json")
  strong_pickle = os.path.join(
    dials_regression, "indexing_test_data/phi_scan/strong.pickle")

  from dxtbx.serialize import load
  imageset_old = load.datablock(
    datablock_old, check_format=False)[0].extract_imagesets()[0]
  imageset_new = load.datablock(
    datablock_new, check_format=False)[0].extract_imagesets()[0]

  gonio_old = imageset_old.get_goniometer()
  gonio_new = imageset_new.get_goniometer()

  assert approx_equal(
    gonio_old.get_rotation_axis(),
    (0.7497646259807715, -0.5517923303436749, 0.36520984351713554))
  assert approx_equal(
    gonio_old.get_setting_rotation(),
    (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
  assert approx_equal(
    gonio_old.get_fixed_rotation(),
    (0.7497646259807748, -0.20997265900532208, -0.6275065641872948,
     -0.5517923303436731, 0.3250014637526764, -0.7680490041218182,
     0.3652098435171313, 0.9221092836691605, 0.12781329809272568))

  assert approx_equal(
    gonio_new.get_rotation_axis(), gonio_old.get_rotation_axis())
  assert approx_equal(gonio_new.get_rotation_axis_datum(), (1,0,0))
  assert approx_equal(
    gonio_new.get_setting_rotation(),
    (0.7497646259807705, -0.20997265900532142, -0.6275065641873,
     -0.5517923303436786, 0.3250014637526763, -0.768049004121814,
     0.3652098435171315, 0.9221092836691607, 0.12781329809272335))
  assert approx_equal(
    gonio_new.get_fixed_rotation(),
    (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

  result_old = run_one_indexing(
    pickle_path=strong_pickle, sweep_path=datablock_old,
    extra_args=[],
    expected_unit_cell=expected_unit_cell,
    expected_rmsds=expected_rmsds,
    expected_hall_symbol=' P 1',
    )

  result_new = run_one_indexing(
    pickle_path=strong_pickle, sweep_path=datablock_new,
    extra_args=[],
    expected_unit_cell=expected_unit_cell,
    expected_rmsds=expected_rmsds,
    expected_hall_symbol=' P 1',
    )

  assert approx_equal(result_old.rmsds, result_new.rmsds)
  assert approx_equal(result_old.crystal_model.get_unit_cell().parameters(),
                      result_new.crystal_model.get_unit_cell().parameters())


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
