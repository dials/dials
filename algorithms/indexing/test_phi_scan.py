from __future__ import absolute_import, division, print_function

import os
import pytest
from cctbx import uctbx
from dials.algorithms.indexing.test_index import run_indexing


def test_run(dials_regression, tmpdir):
    expected_unit_cell = uctbx.unit_cell(
        (11.624, 13.550, 30.103, 89.964, 93.721, 90.132)
    )
    expected_rmsds = (0.039, 0.035, 0.002)

    experiments_old = os.path.join(
        dials_regression, "indexing_test_data", "phi_scan", "datablock_old.json"
    )
    experiments_new = os.path.join(
        dials_regression, "indexing_test_data", "phi_scan", "datablock.json"
    )
    strong_pickle = os.path.join(
        dials_regression, "indexing_test_data", "phi_scan", "strong.pickle"
    )

    from dxtbx.serialize import load

    imageset_old = load.experiment_list(
        experiments_old, check_format=False
    ).imagesets()[0]
    imageset_new = load.experiment_list(
        experiments_new, check_format=False
    ).imagesets()[0]

    gonio_old = imageset_old.get_goniometer()
    gonio_new = imageset_new.get_goniometer()

    assert gonio_old.get_rotation_axis() == pytest.approx(
        (0.7497646259807715, -0.5517923303436749, 0.36520984351713554)
    )
    assert gonio_old.get_setting_rotation() == pytest.approx(
        (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    )
    assert gonio_old.get_fixed_rotation() == pytest.approx(
        (
            0.7497646259807748,
            -0.20997265900532208,
            -0.6275065641872948,
            -0.5517923303436731,
            0.3250014637526764,
            -0.7680490041218182,
            0.3652098435171313,
            0.9221092836691605,
            0.12781329809272568,
        )
    )

    assert gonio_new.get_rotation_axis() == pytest.approx(gonio_old.get_rotation_axis())
    assert gonio_new.get_rotation_axis_datum() == pytest.approx((1, 0, 0))
    assert gonio_new.get_setting_rotation() == pytest.approx(
        (
            0.7497646259807705,
            -0.20997265900532142,
            -0.6275065641873,
            -0.5517923303436786,
            0.3250014637526763,
            -0.768049004121814,
            0.3652098435171315,
            0.9221092836691607,
            0.12781329809272335,
        )
    )
    assert gonio_new.get_fixed_rotation() == pytest.approx(
        (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    )

    result_old = run_indexing(
        strong_pickle,
        experiments_old,
        tmpdir,
        extra_args=[],
        expected_unit_cell=expected_unit_cell,
        expected_rmsds=expected_rmsds,
        expected_hall_symbol=" P 1",
    )

    result_new = run_indexing(
        strong_pickle,
        experiments_new,
        tmpdir,
        extra_args=[],
        expected_unit_cell=expected_unit_cell,
        expected_rmsds=expected_rmsds,
        expected_hall_symbol=" P 1",
    )

    assert result_old.rmsds == pytest.approx(result_new.rmsds, abs=1e-6)
    assert result_old.experiments[
        0
    ].crystal.get_unit_cell().parameters() == pytest.approx(
        result_new.experiments[0].crystal.get_unit_cell().parameters(), abs=1e-6
    )

    # Now test refinement gradients are correct
    from dxtbx.model.experiment_list import ExperimentList, Experiment

    old_exps = ExperimentList(
        [
            Experiment(
                beam=imageset_old.get_beam(),
                detector=imageset_old.get_detector(),
                goniometer=gonio_old,
                scan=imageset_old.get_scan(),
                crystal=result_old.experiments[0].crystal,
                imageset=None,
            )
        ]
    )
    new_exps = ExperimentList(
        [
            Experiment(
                beam=imageset_new.get_beam(),
                detector=imageset_new.get_detector(),
                goniometer=gonio_new,
                scan=imageset_new.get_scan(),
                crystal=result_new.experiments[0].crystal,
                imageset=None,
            )
        ]
    )

    from libtbx.phil import parse
    from dials.algorithms.refinement.refiner import phil_scope

    params = phil_scope.fetch(source=parse("")).extract()
    from dials.algorithms.refinement.refiner import RefinerFactory

    refiner_old = RefinerFactory.from_parameters_data_experiments(
        params, result_old.indexed_reflections, old_exps, verbosity=0
    )
    refiner_new = RefinerFactory.from_parameters_data_experiments(
        params, result_new.indexed_reflections, new_exps, verbosity=0
    )

    # Analytical gradients should be approximately the same in either case
    an_grads_old = refiner_old._pred_param.get_gradients(refiner_old.get_matches())
    an_grads_new = refiner_new._pred_param.get_gradients(refiner_new.get_matches())
    for g1, g2 in zip(an_grads_old, an_grads_new):
        assert g1["dX_dp"] == pytest.approx(g2["dX_dp"], abs=1.0e-6)
        assert g1["dY_dp"] == pytest.approx(g2["dY_dp"], abs=1.0e-6)
        assert g1["dphi_dp"] == pytest.approx(g2["dphi_dp"], abs=1.0e-6)

    # Analytical gradients should be approximately equal to finite difference
    # gradients in either case
    fd_grads_old = calc_fd_grads(refiner_old)
    for g1, g2 in zip(fd_grads_old, an_grads_old):
        assert g1["dX_dp"] == pytest.approx(g2["dX_dp"], abs=5.0e-6)
        assert g1["dY_dp"] == pytest.approx(g2["dY_dp"], abs=5.0e-6)
        assert g1["dphi_dp"] == pytest.approx(g2["dphi_dp"], abs=5.0e-6)
    fd_grads_new = calc_fd_grads(refiner_new)
    for g1, g2 in zip(fd_grads_new, an_grads_new):
        assert g1["dX_dp"] == pytest.approx(g2["dX_dp"], abs=5.0e-6)
        assert g1["dY_dp"] == pytest.approx(g2["dY_dp"], abs=5.0e-6)
        assert g1["dphi_dp"] == pytest.approx(g2["dphi_dp"], abs=5.0e-6)


def calc_fd_grads(refiner):

    p_vals = refiner._pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    fd_grads = []
    for i, val in enumerate(p_vals):
        p_vals[i] -= deltas[i] / 2.0
        refiner._pred_param.set_param_vals(p_vals)

        refiner._target.predict()

        rev_state = refiner.get_matches()["xyzcal.mm"].deep_copy()

        p_vals[i] += deltas[i]
        refiner._pred_param.set_param_vals(p_vals)

        refiner._target.predict()

        fwd_state = refiner.get_matches()["xyzcal.mm"].deep_copy()
        p_vals[i] = val

        fd = fwd_state - rev_state
        x_grads, y_grads, phi_grads = fd.parts()
        x_grads /= deltas[i]
        y_grads /= deltas[i]
        phi_grads /= deltas[i]

        fd_grads.append({"dX_dp": x_grads, "dY_dp": y_grads, "dphi_dp": phi_grads})

    # return to the initial state
    refiner._pred_param.set_param_vals(p_vals)

    return fd_grads
