from __future__ import absolute_import, division, print_function

import pytest
from dials.algorithms.refinement.reflection_manager import phil_scope as refman_phil_scope
from dials.algorithms.refinement.reflection_manager import ReflectionManagerFactory
from dials.algorithms.refinement.parameterisation.autoreduce import phil_scope as ar_phil_scope
from dials.algorithms.refinement.parameterisation.autoreduce import AutoReduce

from dials.test.algorithms.refinement.test_stills_prediction_parameters import _Test
from dials.algorithms.refinement.prediction.managed_predictors import StillsExperimentsPredictor

@pytest.fixture(scope="session")
def tc():
  test = _Test()

  # Predict the reflections in place and put in a reflection manager
  ref_predictor = StillsExperimentsPredictor(test.stills_experiments)
  ref_predictor(test.reflections)
  test.refman = ReflectionManagerFactory.from_parameters_reflections_experiments(
      refman_phil_scope.extract(), test.reflections, test.stills_experiments,
      do_stills=True)
  test.refman.finalise()
  return test

def test_check_and_fail(tc):

  # There are 823 reflections and the detector parameterisation has 6 free
  # parameters
  assert len(tc.refman.get_matches()) == 823
  assert tc.det_param.num_free() == 6

  # Setting 137 reflections as the minimum should pass (137*6<823)
  options = ar_phil_scope.extract()
  options.min_nref_per_parameter=137
  ar = AutoReduce(options, [tc.det_param], [tc.s0_param],
      [tc.xlo_param], [tc.xluc_param], gon_params=[],
      reflection_manager=tc.refman)

  ar.check_and_fail()

  # Setting 138 reflections as the minimum should fail (138*6>823)
  options.min_nref_per_parameter=138
  ar = AutoReduce(options, [tc.det_param], [tc.s0_param],
      [tc.xlo_param], [tc.xluc_param], gon_params=[],
      reflection_manager=tc.refman)
  from libtbx.test_utils import Sorry
  try:
    ar.check_and_fail()
  except Sorry:
    return
  raise RuntimeError

def test_check_and_fix(tc):

  n_det = tc.det_param.num_free()
  n_beam = tc.s0_param.num_free()
  n_xlo = tc.xlo_param.num_free()
  n_xluc = tc.xluc_param.num_free()

  # Similar to test_check_and_fail, setting 137 reflections as the minimum
  # should leave all parameters free
  options = ar_phil_scope.extract()
  options.min_nref_per_parameter=137
  ar = AutoReduce(options, [tc.det_param], [tc.s0_param],
      [tc.xlo_param], [tc.xluc_param], gon_params=[],
      reflection_manager=tc.refman)
  ar.check_and_fix()

  assert ar.det_params[0].num_free() == n_det == 6
  assert ar.beam_params[0].num_free() == n_beam == 3
  assert ar.xl_ori_params[0].num_free() == n_xlo == 3
  assert ar.xl_uc_params[0].num_free() == n_xluc == 6

  # Setting 138 reflections as the minimum should fix all the detector
  # parameters and remove that parameterisation. The crystal unit cell also
  # has 6 parameters, but each parameter is considered separately, so the
  # critical minimum number of reflections is 138*1 not 138*6 in that case
  options = ar_phil_scope.extract()
  options.min_nref_per_parameter=138
  ar = AutoReduce(options, [tc.det_param], [tc.s0_param],
      [tc.xlo_param], [tc.xluc_param], gon_params=[],
      reflection_manager=tc.refman)
  ar.check_and_fix()

  assert not ar.det_params
  assert ar.xl_uc_params[0].num_free() == n_xluc
  assert ar.beam_params[0].num_free() == n_beam
  assert ar.xl_ori_params[0].num_free() == n_xlo

def test_check_and_remove():

  test = _Test()

  # Override the single panel model and parameterisation. This test function
  # exercises the code for non-hierarchical multi-panel detectors. The
  # hierarchical detector version is tested via test_cspad_refinement.py
  from dxtbx.model import Detector
  from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationMultiPanel
  from dials.test.algorithms.refinement.test_multi_panel_detector_parameterisation import \
    make_panel_in_array
  multi_panel_detector = Detector()
  for x in range(3):
    for y in range(3):
      new_panel = make_panel_in_array((x, y), test.detector[0])
      multi_panel_detector.add_panel(new_panel)
  test.detector = multi_panel_detector
  test.stills_experiments[0].detector = multi_panel_detector
  test.det_param = DetectorParameterisationMultiPanel(multi_panel_detector, test.beam)

  # update the generated reflections
  test.generate_reflections()

  # Predict the reflections in place and put in a reflection manager
  ref_predictor = StillsExperimentsPredictor(test.stills_experiments)
  ref_predictor(test.reflections)
  test.refman = ReflectionManagerFactory.from_parameters_reflections_experiments(
      refman_phil_scope.extract(), test.reflections, test.stills_experiments,
      do_stills=True)
  test.refman.finalise()

  # A non-hierarchical detector does not have panel groups, thus panels are
  # not treated independently wrt which reflections affect their parameters.
  # As before, setting 137 reflections as the minimum should leave all
  # parameters free, and should not remove any reflections
  options = ar_phil_scope.extract()
  options.min_nref_per_parameter=137
  ar = AutoReduce(options, [test.det_param], [test.s0_param],
      [test.xlo_param], [test.xluc_param], gon_params=[],
      reflection_manager=test.refman)
  ar.check_and_remove()

  assert ar.det_params[0].num_free() ==  6
  assert ar.beam_params[0].num_free() == 3
  assert ar.xl_ori_params[0].num_free() == 3
  assert ar.xl_uc_params[0].num_free() == 6
  assert len(ar.reflection_manager.get_obs()) == 823

  # Setting reflections as the minimum should fix the detector parameters,
  # which removes that parameterisation. Because all reflections are recorded
  # on that detector, they will all be removed as well. This then affects all
  # other parameterisations, which will be removed.
  options = ar_phil_scope.extract()
  options.min_nref_per_parameter=138
  ar = AutoReduce(options, [test.det_param], [test.s0_param],
      [test.xlo_param], [test.xluc_param], gon_params=[],
      reflection_manager=test.refman)
  ar.check_and_remove()

  assert not ar.det_params
  assert not ar.beam_params
  assert not ar.xl_ori_params
  assert not ar.xl_uc_params
  assert len(ar.reflection_manager.get_obs()) == 0

