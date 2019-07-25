from __future__ import absolute_import, division, print_function

import copy

import pytest
from dials.algorithms.refinement.reflection_manager import (
    phil_scope as refman_phil_scope,
)
from dials.algorithms.refinement.reflection_manager import ReflectionManagerFactory
from dials.algorithms.refinement.parameterisation.autoreduce import (
    phil_scope as ar_phil_scope,
)
from dials.algorithms.refinement.parameterisation.autoreduce import AutoReduce

from dials.test.algorithms.refinement.test_stills_prediction_parameters import _Test
from dials.algorithms.refinement.prediction.managed_predictors import (
    StillsExperimentsPredictor,
)
from dials.array_family import flex
from dials.algorithms.refinement import DialsRefineConfigError


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
    return test


def test_check_and_fail(tc):

    # There are 823 reflections and the detector parameterisation has 6 free
    # parameters
    assert len(tc.refman.get_matches()) == 823
    assert tc.det_param.num_free() == 6

    # Setting 137 reflections as the minimum should pass (137*6<823)
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 137
    ar = AutoReduce(
        options,
        [tc.det_param],
        [tc.s0_param],
        [tc.xlo_param],
        [tc.xluc_param],
        gon_params=[],
        reflection_manager=tc.refman,
    )

    ar.check_and_fail()

    # Setting 138 reflections as the minimum should fail (138*6>823)
    options.min_nref_per_parameter = 138
    ar = AutoReduce(
        options,
        [tc.det_param],
        [tc.s0_param],
        [tc.xlo_param],
        [tc.xluc_param],
        gon_params=[],
        reflection_manager=tc.refman,
    )

    with pytest.raises(DialsRefineConfigError):
        ar.check_and_fail()


def test_check_and_fix(tc):

    n_det = tc.det_param.num_free()
    n_beam = tc.s0_param.num_free()
    n_xlo = tc.xlo_param.num_free()
    n_xluc = tc.xluc_param.num_free()

    # Similar to test_check_and_fail, setting 137 reflections as the minimum
    # should leave all parameters free
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 137
    ar = AutoReduce(
        options,
        [tc.det_param],
        [tc.s0_param],
        [tc.xlo_param],
        [tc.xluc_param],
        gon_params=[],
        reflection_manager=tc.refman,
    )
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
    options.min_nref_per_parameter = 138
    ar = AutoReduce(
        options,
        [tc.det_param],
        [tc.s0_param],
        [tc.xlo_param],
        [tc.xluc_param],
        gon_params=[],
        reflection_manager=tc.refman,
    )
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
    from dials.algorithms.refinement.parameterisation.detector_parameters import (
        DetectorParameterisationMultiPanel,
    )
    from dials.test.algorithms.refinement.test_multi_panel_detector_parameterisation import (
        make_panel_in_array,
    )

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

    # A non-hierarchical detector does not have panel groups, thus panels are
    # not treated independently wrt which reflections affect their parameters.
    # As before, setting 137 reflections as the minimum should leave all
    # parameters free, and should not remove any reflections
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 137
    ar = AutoReduce(
        options,
        [test.det_param],
        [test.s0_param],
        [test.xlo_param],
        [test.xluc_param],
        gon_params=[],
        reflection_manager=test.refman,
    )
    ar.check_and_remove()

    assert ar.det_params[0].num_free() == 6
    assert ar.beam_params[0].num_free() == 3
    assert ar.xl_ori_params[0].num_free() == 3
    assert ar.xl_uc_params[0].num_free() == 6
    assert len(ar.reflection_manager.get_obs()) == 823

    # Setting reflections as the minimum should fix the detector parameters,
    # which removes that parameterisation. Because all reflections are recorded
    # on that detector, they will all be removed as well. This then affects all
    # other parameterisations, which will be removed.
    options = ar_phil_scope.extract()
    options.min_nref_per_parameter = 138
    ar = AutoReduce(
        options,
        [test.det_param],
        [test.s0_param],
        [test.xlo_param],
        [test.xluc_param],
        gon_params=[],
        reflection_manager=test.refman,
    )
    ar.check_and_remove()

    assert not ar.det_params
    assert not ar.beam_params
    assert not ar.xl_ori_params
    assert not ar.xl_uc_params
    assert len(ar.reflection_manager.get_obs()) == 0


# Test the functionality of the parameter 'auto reduction' extension modules
@pytest.fixture(scope="session")
def setup_test_sorting():
    # Borrowed from tst_reflection_table function tst_find_overlapping

    N = 110
    r = flex.reflection_table.empty_standard(N)
    r["panel"] = flex.size_t([1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0] * 10)
    r["id"] = flex.int([1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 1] * 10)
    exp_ids = flex.size_t([0, 1])
    for i in range(N):
        r["miller_index"][i] = (
            int(i // 10) - 5,
            i % 3,
            i % 7,
        )  # A nice bunch of miller indices

    # Filter out reflections to be used by refinement. Sorting of filtered reflections
    # require to allow C++ extension modules to give performance benefit. Sorting
    # performed within the _filter_reflections step by id, then by panel.
    r_sorted = copy.deepcopy(r)
    r_sorted.sort("id")
    r_sorted.subsort("id", "panel")

    # Test that the unfiltered/unsorted table becomes filtered/sorted for id
    assert (r_sorted["id"] == r["id"].select(flex.sort_permutation(r["id"]))).count(
        False
    ) == 0
    # as above for panel within each id
    for ii in [0, 1, 2]:
        r_id = r.select(r["id"] == ii)
        r_sorted_id = r_sorted.select(r_sorted["id"] == ii)
        assert (
            r_sorted_id["panel"]
            == r_id["panel"].select(flex.sort_permutation(r_id["panel"]))
        ).count(False) == 0
    return (r, r_sorted, exp_ids)


def test_auto_reduction_parameter_extension_modules_part1(setup_test_sorting):
    # Cut-down original algorithm for AutoReduce._surplus_reflections

    from dials_refinement_helpers_ext import surpl_iter as surpl

    r, r_sorted, exp_ids = setup_test_sorting
    isel = flex.size_t()
    for exp_id in exp_ids:
        isel.extend((r["id"] == exp_id).iselection())
    res0 = len(isel)

    # Updated algorithm for _surplus_reflections, with templated id column for int and size_t
    res1_unsrt_int = surpl(r["id"], exp_ids).result
    res1_int = surpl(r_sorted["id"], exp_ids).result
    res1_sizet = surpl(flex.size_t(list(r_sorted["id"])), exp_ids).result

    # Check that unsorted list fails, while sorted succeeds for both int and size_t array types
    assert res0 != res1_unsrt_int
    assert res0 == res1_int
    assert res0 == res1_sizet


def test_auto_reduction_parameter_extension_modules_part2(setup_test_sorting):
    # Cut-down original algorithm for AutoReduce._unit_cell_surplus_reflections

    from dials_refinement_helpers_ext import uc_surpl_iter as uc_surpl

    r, r_sorted, exp_ids = setup_test_sorting
    isel = flex.size_t()
    for exp_id in exp_ids:
        isel.extend((r["id"] == exp_id).iselection())
    ref = r.select(isel)
    h = ref["miller_index"].as_vec3_double()
    dB_dp = flex.mat3_double([(1, 2, 3, 4, 5, 6, 7, 8, 9), (0, 1, 0, 1, 0, 1, 0, 1, 0)])
    nref_each_param = []
    for der in dB_dp:
        tst = (der * h).norms()
        nref_each_param.append((tst > 0.0).count(True))
    res0 = min(nref_each_param)

    # Updated algorithm for _unit_cell_surplus_reflections
    res1_unsrt_int = uc_surpl(r["id"], r["miller_index"], exp_ids, dB_dp).result
    res1_int = uc_surpl(r_sorted["id"], r_sorted["miller_index"], exp_ids, dB_dp).result
    res1_sizet = uc_surpl(
        flex.size_t(list(r_sorted["id"])), r_sorted["miller_index"], exp_ids, dB_dp
    ).result
    assert res0 != res1_unsrt_int
    assert res0 == res1_int
    assert res0 == res1_sizet


def test_auto_reduction_parameter_extension_modules_part3(setup_test_sorting):
    # Cut-down original algorithm for AutoReduce._panel_gp_surplus_reflections

    from dials_refinement_helpers_ext import pg_surpl_iter as pg_surpl

    r, r_sorted, exp_ids = setup_test_sorting
    isel = flex.size_t()
    pnl_ids = [0, 1]
    for exp_id in exp_ids:
        sub_expID = (r["id"] == exp_id).iselection()
        sub_panels_expID = r["panel"].select(sub_expID)
        for pnl in pnl_ids:
            isel.extend(sub_expID.select(sub_panels_expID == pnl))
    res0 = len(isel)

    # Updated algorithm for _panel_gp_surplus_reflections
    res1_unsrt_int = pg_surpl(r["id"], r["panel"], pnl_ids, exp_ids, 0).result
    res1_int = pg_surpl(r_sorted["id"], r_sorted["panel"], pnl_ids, exp_ids, 0).result
    res1_sizet = pg_surpl(
        flex.size_t(list(r_sorted["id"])), r_sorted["panel"], pnl_ids, exp_ids, 0
    ).result
    assert res0 != res1_unsrt_int
    assert res0 == res1_int
    assert res0 == res1_sizet
