from __future__ import annotations

import pytest

from dxtbx.model import Detector

from dials.algorithms.refinement import DialsRefineConfigError
from dials.algorithms.refinement.parameterisation.autoreduce import AutoReduce
from dials.algorithms.refinement.parameterisation.autoreduce import (
    phil_scope as ar_phil_scope,
)
from dials.algorithms.refinement.parameterisation.detector_parameters import (
    DetectorParameterisationMultiPanel,
)
from dials.algorithms.refinement.parameterisation.prediction_parameters_stills import (
    StillsPredictionParameterisation,
)
from dials.algorithms.refinement.prediction.managed_predictors import (
    StillsExperimentsPredictor,
)
from dials.algorithms.refinement.reflection_manager import ReflectionManagerFactory
from dials.algorithms.refinement.reflection_manager import (
    phil_scope as refman_phil_scope,
)

from .test_multi_panel_detector_parameterisation import make_panel_in_array
from .test_stills_prediction_parameters import _Test


@pytest.fixture(scope="session")
def tc():
    test = _Test()

    # Predict the reflections in place and put in a reflection manager
    ref_predictor = StillsExperimentsPredictor(test.stills_experiments)
    ref_predictor(test.reflections)
    test.refman = ReflectionManagerFactory.from_parameters_reflections_experiments(
        refman_phil_scope.extract(),
        test.reflections,
        test.stills_experiments,
        do_stills=True,
    )
    test.refman.finalise()

    # Build a prediction parameterisation for the stills experiment
    test.pred_param = StillsPredictionParameterisation(
        test.stills_experiments,
        detector_parameterisations=[test.det_param],
        beam_parameterisations=[test.s0_param],
        xl_orientation_parameterisations=[test.xlo_param],
        xl_unit_cell_parameterisations=[test.xluc_param],
    )
    return test


def test_check_and_fail(tc):

    # There are 823 reflections
    assert len(tc.refman.get_matches()) == 823

    # The parameters affecting the smallest number of reflections are
    # g_param_0, g_param_3 and g_param_4, all of which have gradients for
    # 792 reflections. Setting 792 reflections as the minimum should pass.
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 792
    ar = AutoReduce(options, pred_param=tc.pred_param, reflection_manager=tc.refman)

    ar.check_and_fail()

    # Setting 793 reflections as the minimum should fail
    options.min_nref_per_parameter = 793
    ar = AutoReduce(options, pred_param=tc.pred_param, reflection_manager=tc.refman)

    with pytest.raises(DialsRefineConfigError):
        ar.check_and_fail()


def test_check_and_fix(tc):

    n_det = tc.det_param.num_free()
    n_beam = tc.s0_param.num_free()
    n_xlo = tc.xlo_param.num_free()
    n_xluc = tc.xluc_param.num_free()

    # Similar to test_check_and_fail, setting 792 reflections as the minimum
    # should leave all parameters free
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 792
    ar = AutoReduce(options, pred_param=tc.pred_param, reflection_manager=tc.refman)
    ar.check_and_fix()

    det_params = tc.pred_param.get_detector_parameterisations()
    beam_params = tc.pred_param.get_beam_parameterisations()
    xl_ori_params = tc.pred_param.get_crystal_orientation_parameterisations()
    xl_uc_params = tc.pred_param.get_crystal_unit_cell_parameterisations()
    assert det_params[0].num_free() == n_det == 6
    assert beam_params[0].num_free() == n_beam == 3
    assert xl_ori_params[0].num_free() == n_xlo == 3
    assert xl_uc_params[0].num_free() == n_xluc == 6

    # Setting 793 reflections as the minimum should fix crystal unit cell
    # parameters g_param_0, g_param_3 and g_param_4
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 793
    ar = AutoReduce(options, pred_param=tc.pred_param, reflection_manager=tc.refman)
    ar.check_and_fix()

    det_params = tc.pred_param.get_detector_parameterisations()
    beam_params = tc.pred_param.get_beam_parameterisations()
    xl_ori_params = tc.pred_param.get_crystal_orientation_parameterisations()
    xl_uc_params = tc.pred_param.get_crystal_unit_cell_parameterisations()
    assert det_params[0].num_free() == n_det
    assert xl_uc_params[0].num_free() == 3
    assert xl_uc_params[0].get_param_names() == ["g_param_1", "g_param_2", "g_param_5"]
    assert beam_params[0].num_free() == n_beam
    assert xl_ori_params[0].num_free() == n_xlo


def test_check_and_remove():

    test = _Test()

    # Override the single panel model and parameterisation. This test function
    # exercises the code for non-hierarchical multi-panel detectors. The
    # hierarchical detector version is tested via test_cspad_refinement.py
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
        refman_phil_scope.extract(),
        test.reflections,
        test.stills_experiments,
        do_stills=True,
    )
    test.refman.finalise()

    # Build a prediction parameterisation for the stills experiment
    test.pred_param = StillsPredictionParameterisation(
        test.stills_experiments,
        detector_parameterisations=[test.det_param],
        beam_parameterisations=[test.s0_param],
        xl_orientation_parameterisations=[test.xlo_param],
        xl_unit_cell_parameterisations=[test.xluc_param],
    )

    # A non-hierarchical detector does not have panel groups, thus panels are
    # not treated independently wrt which reflections affect their parameters.
    # As before, setting 792 reflections as the minimum should leave all
    # parameters free, and should not remove any reflections
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 792
    ar = AutoReduce(options, pred_param=test.pred_param, reflection_manager=test.refman)
    ar.check_and_remove()

    det_params = test.pred_param.get_detector_parameterisations()
    beam_params = test.pred_param.get_beam_parameterisations()
    xl_ori_params = test.pred_param.get_crystal_orientation_parameterisations()
    xl_uc_params = test.pred_param.get_crystal_unit_cell_parameterisations()
    assert det_params[0].num_free() == 6
    assert beam_params[0].num_free() == 3
    assert xl_ori_params[0].num_free() == 3
    assert xl_uc_params[0].num_free() == 6
    assert len(test.refman.get_obs()) == 823

    # Setting 793 reflections as the minimum fixes 3 unit cell parameters,
    # and removes all those reflections. There are then too few reflections
    # for any parameterisation and all will be fixed, leaving no free
    # parameters for refinement. This fails within PredictionParameterisation,
    # during update so the final 31 reflections are not removed.
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 793
    ar = AutoReduce(options, pred_param=test.pred_param, reflection_manager=test.refman)
    with pytest.raises(
        DialsRefineConfigError, match="There are no free parameters for refinement"
    ):
        ar.check_and_remove()

    det_params = test.pred_param.get_detector_parameterisations()
    beam_params = test.pred_param.get_beam_parameterisations()
    xl_ori_params = test.pred_param.get_crystal_orientation_parameterisations()
    xl_uc_params = test.pred_param.get_crystal_unit_cell_parameterisations()
    assert det_params[0].num_free() == 0
    assert beam_params[0].num_free() == 0
    assert xl_ori_params[0].num_free() == 0
    assert xl_uc_params[0].num_free() == 0
    assert len(test.refman.get_obs()) == 823 - 792
