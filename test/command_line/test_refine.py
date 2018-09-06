"""
Test command line program dials.refine by running a job with saved data and
comparing with expected output.

This serves as a high level test that not only checks whether refinement works,
but also that the command line program is functioning and that the output models
have not changed format and so on.
"""

from __future__ import absolute_import, division, print_function

import copy
import six.moves.cPickle as pickle
import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
import pytest

def test1(dials_regression, tmpdir):
  tmpdir.chdir()

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  experiments_path = os.path.join(data_dir, "experiments.json")
  pickle_path = os.path.join(data_dir, "indexed_strong.pickle")

  for pth in (experiments_path, pickle_path):
    assert os.path.exists(pth)

  # set some old defaults
  cmd = "dials.refine close_to_spindle_cutoff=0.05 reflections_per_degree=100 " + \
        "outlier.separate_blocks=False " + experiments_path + " " + pickle_path

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  # load results
  reg_exp = ExperimentListFactory.from_json_file(
              os.path.join(data_dir, "regression_experiments.json"),
              check_format=False)[0]
  ref_exp = ExperimentListFactory.from_json_file("refined_experiments.json",
              check_format=False)[0]

  # test refined models against expected
  assert reg_exp.crystal == ref_exp.crystal
  assert reg_exp.detector == ref_exp.detector
  assert reg_exp.beam == ref_exp.beam

  # test cell parameter esds
  assert ref_exp.crystal.get_cell_parameter_sd() == \
      pytest.approx((0.0009903, 0.0009903, 0.0021227, 0.0, 0.0, 0.0), abs=1e-6)
  assert ref_exp.crystal.get_cell_volume_sd() == pytest.approx(23.8063382, abs=1e-6)

def test2(dials_regression, tmpdir):
  """Run scan-varying refinement, comparing RMSD table with expected values.
  This test automates what was manually done periodically and recorded in
  dials_regression/refinement_test_data/centroid/README.txt"""

  tmpdir.chdir()

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "centroid")
  experiments_path = os.path.join(data_dir, "experiments_XPARM_REGULARIZED.json")
  pickle_path = os.path.join(data_dir, "spot_all_xds.pickle")

  for pth in (experiments_path, pickle_path):
    assert os.path.exists(pth)

  # scan-static refinement first to get refined_experiments.json as start point
  cmd1 = ("dials.refine " + experiments_path + " " + pickle_path +
    " reflections_per_degree=50 "
    " outlier.algorithm=null close_to_spindle_cutoff=0.05")
  cmd2 = ("dials.refine refined_experiments.json " + pickle_path +
    " scan_varying=true output.history=history.pickle"
    " reflections_per_degree=50"
    " outlier.algorithm=null close_to_spindle_cutoff=0.05"
    " crystal.orientation.smoother.interval_width_degrees=36.0"
    " crystal.unit_cell.smoother.interval_width_degrees=36.0"
    " set_scan_varying_errors=True")

  result1 = easy_run.fully_buffered(command=cmd1).raise_if_errors()
  result2 = easy_run.fully_buffered(command=cmd2).raise_if_errors()

  # load and check results
  with open("history.pickle", "rb") as fh:
    history=pickle.load(fh)

  expected_rmsds = [(0.088488398, 0.114583571, 0.001460382),
                    (0.080489334, 0.086406517, 0.001284069),
                    (0.078835086, 0.086052630, 0.001195882),
                    (0.077476911, 0.086194611, 0.001161143),
                    (0.076755840, 0.086090630, 0.001157239),
                    (0.076586376, 0.085939462, 0.001155641),
                    (0.076603722, 0.085878953, 0.001155065),
                    (0.076611382, 0.085862959, 0.001154863),
                    (0.076608732, 0.085856935, 0.001154384),
                    (0.076605731, 0.085852271, 0.001153858),
                    (0.076604576, 0.085852318, 0.001153643),
                    (0.076603981, 0.085854175, 0.001153594)]
  assert approx_equal(history['rmsd'], expected_rmsds)

  # check that the used_in_refinement flag got set correctly
  rt = flex.reflection_table.from_pickle('refined.pickle')
  uir = rt.get_flags(rt.flags.used_in_refinement)
  assert uir.count(True) == history['num_reflections'][-1]

def test3(dials_regression, tmpdir):
  """Strict check for scan-varying refinement using automated outlier rejection
  block width and interval width setting"""
  tmpdir.chdir()

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "centroid")
  experiments_path = os.path.join(data_dir, "experiments_XPARM_REGULARIZED.json")
  pickle_path = os.path.join(data_dir, "spot_all_xds.pickle")

  for pth in (experiments_path, pickle_path):
    assert os.path.exists(pth)

  cmd1 = ("dials.refine " + experiments_path + " " + pickle_path +
          " scan_varying=true max_iterations=5 output.history=history.pickle "
          "crystal.orientation.smoother.interval_width_degrees=auto "
          "crystal.unit_cell.smoother.interval_width_degrees=auto")
  result1 = easy_run.fully_buffered(command=cmd1).raise_if_errors()

  # load and check results
  with open("history.pickle", "rb") as fh:
    history=pickle.load(fh)

  expected_rmsds = [[0.619507829, 0.351326044, 0.006955399],
                    [0.174024575, 0.113486044, 0.004704006],
                    [0.098351363, 0.084052519, 0.002660408],
                    [0.069202909, 0.072796782, 0.001451734],
                    [0.064305277, 0.071560831, 0.001165639],
                    [0.062955462, 0.071315612, 0.001074453]]
  assert approx_equal(history['rmsd'], expected_rmsds)

  # check the refined unit cell
  ref_exp = ExperimentListFactory.from_json_file("refined_experiments.json",
    check_format=False)[0]
  unit_cell = ref_exp.crystal.get_unit_cell().parameters()
  assert unit_cell == pytest.approx([42.27482, 42.27482, 39.66893, 90.00000, 90.00000, 90.00000],
      abs=1e-3)


#Test the functionality of the refiner.py extension modules
def test4():
  from dials_refinement_helpers_ext import pgnmn_iter as pgnmn
  from dials_refinement_helpers_ext import ucnmn_iter as ucnmn
  from dials_refinement_helpers_ext import mnmn_iter as mnmn
  #Borrowed from tst_reflection_table function tst_find_overlapping
  from random import randint, uniform
  N = 110
  r = flex.reflection_table.empty_standard(N)
  r['panel'] = flex.size_t([1,0,0,1,0,0,0,0,1,0,0]*10)
  r['id'] = flex.int([1,2,1,1,2,0,1,1,1,0,1]*10)
  exp_ids = flex.size_t([0,1])
  for i in xrange(N):
    r['miller_index'][i] = (int(i//10) - 5, i%3, i%7) #A nice bunch of miller indices

  '''
   Filter out reflections to be used by refinement. Sorting of filtered reflections require
   to allow C++ extension modules to give performance benefit. Sorting performed within the
   _filter_reflections step by id, then by panel.
  '''
  r_sorted = copy.deepcopy(r)
  r_sorted.sort('id')
  r_sorted.subsort('id','panel')

  # Test that the unfiltered/unsorted table becomes filtered/sorted for id
  assert (r_sorted['id']==r['id'].select(flex.sort_permutation(r['id']))).count(False) == 0
  # as above for panel within each id
  for ii in [0,1,2]:
    r_id = r.select(r['id']==ii)
    r_sorted_id = r_sorted.select(r_sorted['id']==ii)
    assert (r_sorted_id['panel']==r_id['panel'].select(flex.sort_permutation(r_id['panel']))).count(False) == 0

  ############################################################
  #Cut-down original algorithm for model_nparam_minus_nref
  ############################################################
  isel = flex.size_t()
  for exp_id in exp_ids:
    isel.extend((r['id'] == exp_id).iselection())
  res0 = len(isel)

  #Updated algorithm for model_nparam_minus_nref, with templated id column for int and size_t
  res1_unsrt_int = mnmn(r["id"],exp_ids).result
  res1_int = mnmn(r_sorted["id"],exp_ids).result
  res1_sizet = mnmn(flex.size_t(list(r_sorted["id"])),exp_ids).result

  #Check that unsorted list fails, while sorted succeeds for both int and size_t array types
  assert res0 != res1_unsrt_int
  assert res0 == res1_int
  assert res0 == res1_sizet

  ############################################################
  #Cut-down original algorithm for unit_cell_nparam_minus_nref
  ############################################################
  ref = r_sorted.select(isel)
  h = ref['miller_index'].as_vec3_double()
  dB_dp = flex.mat3_double( [(1,2,3,4,5,6,7,8,9), (0,1,0,1,0,1,0,1,0) ] )
  nref_each_param = []
  for der in dB_dp:
    tst = (der * h).norms()
    nref_each_param.append((tst > 0.0).count(True))
  res0 = min(nref_each_param)

  #Updated algorithm for unit_cell_nparam_minus_nref
  res1_unsrt_int = ucnmn(r["id"], r["miller_index"], exp_ids, dB_dp).result
  res1_int = ucnmn(r_sorted["id"], r_sorted["miller_index"], exp_ids, dB_dp).result
  res1_sizet = ucnmn(flex.size_t(list(r_sorted["id"])), r_sorted["miller_index"], exp_ids, dB_dp).result
  assert res0 != res1_unsrt_int
  assert res0 == res1_int
  assert res0 == res1_sizet

  ############################################################
  #Cut-down original algorithm for panel_gp_nparam_minus_nref
  ############################################################
  isel = flex.size_t()
  pnl_ids = [0,1]
  for exp_id in exp_ids:
    sub_expID = (r['id'] == exp_id).iselection()
    sub_panels_expID = r['panel'].select(sub_expID)
    for pnl in pnl_ids:
      isel.extend(sub_expID.select(sub_panels_expID == pnl))
  nref = len(isel)
  res0 = nref

  #Updated algorithm for panel_gp_nparam_minus_nref
  res1_unsrt_int = pgnmn(r["id"], r["panel"], pnl_ids, exp_ids, 0).result
  res1_int = pgnmn(r_sorted["id"], r_sorted["panel"], pnl_ids, exp_ids, 0).result
  res1_sizet = pgnmn(flex.size_t(list(r_sorted["id"])), r_sorted["panel"], pnl_ids, exp_ids, 0).result
  assert res0 != res1_unsrt_int
  assert res0 == res1_int
  assert res0 == res1_sizet
