from __future__ import annotations

from collections import namedtuple
from copy import copy
from random import randint, uniform

import numpy as np
import pytest

from scitbx import matrix

from dials.algorithms.profile_model.ellipsoid.parameterisation import (
    Angular2MosaicityParameterisation,
    Angular4MosaicityParameterisation,
    ModelState,
    ReflectionModelState,
    Simple1MosaicityParameterisation,
    Simple6MosaicityParameterisation,
    WavelengthSpreadParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)
from dials.array_family import flex


def test_Simple1MosaicityParameterisation():

    p = Simple1MosaicityParameterisation(params=np.array([1e-3]))

    assert p.is_angular() is False
    assert p.num_parameters() == 1
    assert p.parameters == (1e-3,)
    p.parameters = np.array([2e-3])
    assert p.parameters == (2e-3,)
    psq = p.parameters[0] ** 2
    assert list(p.sigma().flatten()) == pytest.approx([psq, 0, 0, 0, psq, 0, 0, 0, psq])
    d = p.first_derivatives()
    assert d.shape[0] == 1
    d1 = 2 * p.parameters[0]
    assert list(d[0, :, :].flatten()) == pytest.approx([d1, 0, 0, 0, d1, 0, 0, 0, d1])


def test_Simple6MosaicityParameterisation():

    params = np.array([1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3])

    p = Simple6MosaicityParameterisation(params=params)

    assert p.is_angular() is False
    assert p.num_parameters() == 6
    assert list(p.parameters) == pytest.approx(list(params))

    params = np.array([2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3])
    p.parameters = params
    assert list(p.parameters) == pytest.approx(list(params))

    b1, b2, b3, b4, b5, b6 = params
    assert p.sigma()[0, 0] == pytest.approx(b1**2)
    assert p.sigma()[0, 1] == pytest.approx(b1 * b2)
    assert p.sigma()[0, 2] == pytest.approx(b1 * b4)
    assert p.sigma()[1, 0] == pytest.approx(b1 * b2)
    assert p.sigma()[1, 1] == pytest.approx(b2**2 + b3 * b3)
    assert p.sigma()[1, 2] == pytest.approx(b2 * b4 + b3 * b5)
    assert p.sigma()[2, 0] == pytest.approx(b1 * b4)
    assert p.sigma()[2, 1] == pytest.approx(b2 * b4 + b3 * b5)
    assert p.sigma()[2, 2] == pytest.approx(b4**2 + b5**2 + b6**2)

    dSdb = [
        (2 * b1, b2, b4, b2, 0, 0, b4, 0, 0),
        (0, b1, 0, b1, 2 * b2, b4, 0, b4, 0),
        (0, 0, 0, 0, 2 * b3, b5, 0, b5, 0),
        (0, 0, b1, 0, 0, b2, b1, b2, 2 * b4),
        (0, 0, 0, 0, 0, b3, 0, b3, 2 * b5),
        (0, 0, 0, 0, 0, 0, 0, 0, 2 * b6),
    ]

    d = p.first_derivatives()
    assert d.shape[0] == 6
    for i in range(d.shape[0]):
        a = dSdb[i]
        b = d[i, :, :]
        for j in range(9):
            assert b.flatten()[j] == pytest.approx(a[j], abs=1e-12)


def test_WavelengthSpreadParameterisation():

    params = np.array([1e-3])

    p = WavelengthSpreadParameterisation(params=params)
    assert p.num_parameters() == 1
    assert p.parameters[0] == pytest.approx(params[0])
    params = np.array([2e-3])
    p.parameters = params
    assert p.parameters[0] == pytest.approx(params[0])
    assert p.sigma() == pytest.approx(params[0] ** 2)
    assert p.first_derivatives()[0] == pytest.approx(2 * params[0])


def test_Angular2MosaicityParameterisation():
    params = np.array([1e-3, 2e-3])

    p = Angular2MosaicityParameterisation(params=params)

    assert p.is_angular() is True
    assert p.num_parameters() == 2
    assert p.parameters[0] == pytest.approx(params[0])
    assert p.parameters[1] == pytest.approx(params[1])

    params = np.array([2e-3, 3e-3])
    p.parameters = params
    assert p.parameters[0] == pytest.approx(params[0])
    assert p.parameters[1] == pytest.approx(params[1])

    b1, b2 = params
    assert list(p.sigma().flatten()) == pytest.approx(
        [b1**2, 0, 0, 0, b1**2, 0, 0, 0, b2**2]
    )

    dSdb = [(2 * b1, 0, 0, 0, 2 * b1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 2 * b2)]

    d = p.first_derivatives()
    assert d.shape[0] == 2
    for i in range(d.shape[0]):
        a = dSdb[i]
        b = d[i, :, :]
        for j in range(9):
            assert b.flatten()[j] == pytest.approx(a[j])


def test_Angular4MosaicityParameterisation():
    params = np.array([1e-3, 2e-3, 3e-3, 4e-3])

    p = Angular4MosaicityParameterisation(params=params)

    assert p.is_angular() is True
    assert p.num_parameters() == 4
    assert list(p.parameters) == pytest.approx(list(params))

    params = np.array([2e-3, 3e-3, 4e-3, 5e-3])
    p.parameters = params
    assert list(p.parameters) == pytest.approx(list(params))

    b1, b2, b3, b4 = params
    assert p.sigma()[0, 0] == pytest.approx(b1**2)
    assert p.sigma()[0, 1] == pytest.approx(b1 * b2)
    assert p.sigma()[0, 2] == pytest.approx(0)
    assert p.sigma()[1, 0] == pytest.approx(b1 * b2)
    assert p.sigma()[1, 1] == pytest.approx(b2**2 + b3 * b3)
    assert p.sigma()[1, 2] == pytest.approx(0)
    assert p.sigma()[2, 0] == pytest.approx(0)
    assert p.sigma()[2, 1] == pytest.approx(0)
    assert p.sigma()[2, 2] == pytest.approx(b4**2)

    dSdb = [
        (2 * b1, b2, 0, b2, 0, 0, 0, 0, 0),
        (0, b1, 0, b1, 2 * b2, 0, 0, 0, 0),
        (0, 0, 0, 0, 2 * b3, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 2 * b4),
    ]

    d = p.first_derivatives()
    assert d.shape[0] == 4
    for i in range(d.shape[0]):
        a = dSdb[i]
        b = d[i, ::, :]
        for j in range(9):
            assert b.flatten()[j] == pytest.approx(a[j])


def check_model_state_with_fixed(
    experiment,
    mosaicity_parameterisation,
    wavelength_parameterisation,
    fix_mosaic_spread=False,
    fix_wavelength_spread=False,
    fix_unit_cell=False,
    fix_orientation=False,
):

    state = ModelState(
        experiment,
        mosaicity_parameterisation,
        wavelength_parameterisation,
        fix_mosaic_spread=fix_mosaic_spread,
        fix_wavelength_spread=fix_wavelength_spread,
        fix_unit_cell=fix_unit_cell,
        fix_orientation=fix_orientation,
    )

    assert state.is_orientation_fixed == fix_orientation
    assert state.is_unit_cell_fixed == fix_unit_cell
    assert state.is_mosaic_spread_fixed == fix_mosaic_spread
    assert state.is_wavelength_spread_fixed == fix_wavelength_spread

    U = state.U_matrix
    B = state.B_matrix
    A = state.A_matrix
    M = state.mosaicity_covariance_matrix
    L = state.wavelength_spread

    assert U.shape == (3, 3)
    assert B.shape == (3, 3)
    assert A.shape == (3, 3)
    assert M.shape == (3, 3)
    if wavelength_parameterisation is not None:
        assert len(L) == 1
    else:
        assert len(L) == 0

    assert len(state.U_params) == 3
    assert len(state.B_params) == 2
    assert len(state.M_params) == mosaicity_parameterisation.num_parameters()
    if wavelength_parameterisation is not None:
        assert len(state.L_params) == 1
    else:
        assert len(state.L_params) == 0

    dU = state.dU_dp
    dB = state.dB_dp
    dM = state.dM_dp
    dL = state.dL_dp
    assert dU.shape[0] == 3
    assert dB.shape[0] == 2
    assert dM.shape[0] == mosaicity_parameterisation.num_parameters()
    if wavelength_parameterisation is not None:
        assert len(dL) == 1
    else:
        assert len(dL) == 0

    params = state.active_parameters

    expected_len = 0
    if not fix_mosaic_spread:
        expected_len += mosaicity_parameterisation.num_parameters()
    if not fix_wavelength_spread:
        if wavelength_parameterisation is not None:
            expected_len += wavelength_parameterisation.num_parameters()
    if not fix_unit_cell:
        expected_len += 2
    if not fix_orientation:
        expected_len += 3

    assert len(params) == expected_len
    new_params = params
    state.active_parameters = new_params


def test_ModelState(test_experiment):

    experiments = [test_experiment]

    S1 = Simple1MosaicityParameterisation()
    S6 = Simple6MosaicityParameterisation()
    W = WavelengthSpreadParameterisation()

    with pytest.raises(AssertionError):
        check_model_state_with_fixed(experiments[0], S1, None, fix_mosaic_spread=True)
    with pytest.raises(AssertionError):
        check_model_state_with_fixed(experiments[0], S1, None, fix_unit_cell=True)
    with pytest.raises(AssertionError):
        check_model_state_with_fixed(experiments[0], S1, None, fix_orientation=True)
    check_model_state_with_fixed(experiments[0], S1, None, fix_wavelength_spread=True)
    check_model_state_with_fixed(experiments[0], S1, W, fix_mosaic_spread=True)
    check_model_state_with_fixed(experiments[0], S1, W, fix_wavelength_spread=True)
    check_model_state_with_fixed(experiments[0], S1, W, fix_unit_cell=True)
    check_model_state_with_fixed(experiments[0], S1, W, fix_orientation=True)

    with pytest.raises(AssertionError):
        check_model_state_with_fixed(experiments[0], S6, None, fix_mosaic_spread=True)
    with pytest.raises(AssertionError):
        check_model_state_with_fixed(experiments[0], S6, None, fix_unit_cell=True)
    with pytest.raises(AssertionError):
        check_model_state_with_fixed(experiments[0], S6, None, fix_orientation=True)
    check_model_state_with_fixed(experiments[0], S6, None, fix_wavelength_spread=True)
    check_model_state_with_fixed(experiments[0], S6, W, fix_mosaic_spread=True)
    check_model_state_with_fixed(experiments[0], S6, W, fix_wavelength_spread=True)
    check_model_state_with_fixed(experiments[0], S6, W, fix_unit_cell=True)
    check_model_state_with_fixed(experiments[0], S6, W, fix_orientation=True)


def check_reflection_model_state_with_fixed(
    experiment,
    mosaicity_parameterisation,
    wavelength_parameterisation,
    fix_mosaic_spread=False,
    fix_wavelength_spread=False,
    fix_unit_cell=False,
    fix_orientation=False,
):

    state = ModelState(
        experiment,
        mosaicity_parameterisation,
        wavelength_parameterisation,
        fix_mosaic_spread=fix_mosaic_spread,
        fix_wavelength_spread=fix_wavelength_spread,
        fix_unit_cell=fix_unit_cell,
        fix_orientation=fix_orientation,
    )

    model = ReflectionModelState(
        state, matrix.col(experiment.beam.get_s0()), matrix.col((1, 1, 1))
    )

    assert list(model.mosaicity_covariance_matrix.flatten()) == list(
        mosaicity_parameterisation.sigma().flatten()
    )
    assert list(model.get_r().flatten()) == pytest.approx(
        np.matmul(state.A_matrix, np.array([1, 1, 1]).reshape(3, 1))[:, 0].tolist(),
        abs=1e-6,
    )

    if wavelength_parameterisation is not None:
        assert model.wavelength_spread == wavelength_parameterisation.sigma()
    else:
        assert model.wavelength_spread == 0

    dS_dp = model.get_dS_dp()
    dr_dp = model.get_dr_dp()
    dL_dp = model.get_dL_dp()

    assert dS_dp.shape[2] == len(state.parameter_labels)
    assert dr_dp.shape[1] == len(state.parameter_labels)
    assert len(dL_dp) == len(state.parameter_labels)

    if not fix_wavelength_spread:
        assert dr_dp[:, -1].flatten() == pytest.approx([0, 0, 0], abs=1e-6)
        assert dS_dp[:, :, -1].flatten() == pytest.approx(
            (0, 0, 0, 0, 0, 0, 0, 0, 0), abs=1e-6
        )
        dr_dp = dr_dp[:, :-1]
        dS_dp = dS_dp[:, :, :-1]
        dL_dp = dL_dp[:-1]
    if not fix_mosaic_spread:
        num_params = mosaicity_parameterisation.num_parameters()
        for i in range(num_params):
            assert dr_dp[:, -(i + 1)] == pytest.approx([0, 0, 0], abs=1e-6)
            assert dS_dp[:, :, -(i + 1)] == pytest.approx(
                state.dM_dp[-(i + 1), :, :], abs=1e-6
            )
            assert dL_dp[-1] == 0
        dr_dp = dr_dp[:, :-num_params]
        dS_dp = dS_dp[:, :, :-num_params]
        dL_dp = dL_dp[:-num_params]

    if not fix_unit_cell:
        num_params = state.B_params.size
        for i in range(num_params):
            assert dS_dp[:, :, -(i + 1)].flatten() == pytest.approx(
                (0, 0, 0, 0, 0, 0, 0, 0, 0), abs=1e-6
            )
            assert dL_dp[-(i + 1)] == 0
        dr_dp = dr_dp[:, :-num_params]
        dS_dp = dS_dp[:, :, :-num_params]
        dL_dp = dL_dp[:-num_params]

    if not fix_orientation:
        num_params = state.U_params.size
        for i in range(num_params):
            assert dS_dp[:, :, -(i + 1)].flatten() == pytest.approx(
                (0, 0, 0, 0, 0, 0, 0, 0, 0), abs=1e-6
            )
            assert dL_dp[-(i + 1)] == 0
        dr_dp = dr_dp[:, :-num_params]
        dS_dp = dS_dp[:, :, :-num_params]
        dL_dp = dL_dp[:-num_params]


def test_ReflectionModelState(test_experiment):

    experiments = [test_experiment]

    S1 = Simple1MosaicityParameterisation()
    S6 = Simple6MosaicityParameterisation()
    W = WavelengthSpreadParameterisation()

    check_reflection_model_state_with_fixed(
        experiments[0], S1, None, fix_wavelength_spread=True
    )
    check_reflection_model_state_with_fixed(
        experiments[0], S1, W, fix_mosaic_spread=True
    )
    check_reflection_model_state_with_fixed(
        experiments[0], S1, W, fix_wavelength_spread=True
    )
    check_reflection_model_state_with_fixed(experiments[0], S1, W, fix_unit_cell=True)
    check_reflection_model_state_with_fixed(experiments[0], S1, W, fix_orientation=True)

    check_reflection_model_state_with_fixed(
        experiments[0], S6, None, fix_wavelength_spread=True
    )
    check_reflection_model_state_with_fixed(
        experiments[0], S6, W, fix_mosaic_spread=True
    )
    check_reflection_model_state_with_fixed(
        experiments[0], S6, W, fix_wavelength_spread=True
    )
    check_reflection_model_state_with_fixed(experiments[0], S6, W, fix_unit_cell=True)
    check_reflection_model_state_with_fixed(experiments[0], S6, W, fix_orientation=True)


def first_derivative(func, x, h):
    return (-func(x + 2 * h) + 8 * func(x + h) - 8 * func(x - h) + func(x - 2 * h)) / (
        12 * h
    )


def generate_data(experiments, reflections):

    from random import seed

    seed(0)

    index = randint(0, len(reflections))

    h = reflections[index]["miller_index"]

    s0 = matrix.col(experiments[0].beam.get_s0())

    U_param = CrystalOrientationParameterisation(experiments[0].crystal)
    B_param = CrystalUnitCellParameterisation(experiments[0].crystal)

    U = matrix.sqr(experiments[0].crystal.get_U())
    B = matrix.sqr(experiments[0].crystal.get_B())
    r = U * B * matrix.col(h)
    s2 = s0 + r
    mobs = (
        s2 + matrix.col((uniform(0, 1e-3), uniform(0, 1e-3), uniform(0, 1e-3)))
    ).normalize() * s0.length()

    b1, b2, b3, b4, b5, b6 = (
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
    )

    S_param = (b1, b2, b3, b4, b5, b6)
    L_param = (uniform(1e-3, 2e-3),)
    ctot = randint(100, 1000)

    T = matrix.sqr((uniform(1e-3, 2e-3), 0, uniform(1e-6, 2e-6), uniform(1e-3, 2e-3)))
    Sobs = T * T.transpose()

    params = [S_param, U_param, B_param, L_param]

    return params, s0, h, ctot, mobs, Sobs


@pytest.fixture
def testdata(test_experiment):

    TestData = namedtuple(
        "TestData", ["experiment", "models", "s0", "h", "ctot", "mobs", "Sobs"]
    )

    experiments = [test_experiment]
    reflections = flex.reflection_table.from_predictions_multi(experiments)

    models, s0, h, ctot, mobs, Sobs = generate_data(experiments, reflections)

    return TestData(
        experiment=experiments[0],
        models=models,
        s0=s0,
        h=h,
        ctot=ctot,
        mobs=mobs,
        Sobs=Sobs,
    )


def test_ReflectionModelState_derivatives(testdata):
    def check(
        mosaicity_parameterisation,
        wavelength_parameterisation,
        fix_mosaic_spread=False,
        fix_wavelength_spread=False,
        fix_unit_cell=False,
        fix_orientation=False,
    ):
        experiment = testdata.experiment
        models = testdata.models
        s0 = testdata.s0
        h = testdata.h
        # ctot = testdata.ctot
        # mobs = testdata.mobs
        # Sobs = testdata.Sobs

        U_params = models[1].get_param_vals()
        B_params = models[2].get_param_vals()
        M_params = np.array(models[0][: mosaicity_parameterisation.num_parameters()])
        L_params = np.array(models[3])

        state = ModelState(
            experiment,
            mosaicity_parameterisation,
            wavelength_parameterisation,
            fix_mosaic_spread=fix_mosaic_spread,
            fix_wavelength_spread=fix_wavelength_spread,
            fix_unit_cell=fix_unit_cell,
            fix_orientation=fix_orientation,
        )
        state.U_params = U_params
        state.B_params = B_params
        state.M_params = M_params
        state.L_params = L_params

        model = ReflectionModelState(state, s0, h)

        dr_dp = model.get_dr_dp()
        dS_dp = model.get_dS_dp()
        dL_dp = model.get_dL_dp()

        def compute_sigma(parameters):
            state.active_parameters = parameters
            model = ReflectionModelState(state, s0, h)
            return model.mosaicity_covariance_matrix

        def compute_r(parameters):
            state.active_parameters = parameters
            model = ReflectionModelState(state, s0, h)
            return model.get_r()

        def compute_sigma_lambda(parameters):
            state.active_parameters = parameters
            model = ReflectionModelState(state, s0, h)
            return model.wavelength_spread

        step = 1e-6

        parameters = copy(state.active_parameters)

        dr_num = []
        for i in range(len(parameters)):

            def f(x):
                p = copy(parameters)
                p[i] = x
                return compute_r(p)

            dr_num.append(first_derivative(f, parameters[i], step).reshape(3, 1))
        dr_num = np.concatenate(dr_num, axis=1)

        for n, c in zip(dr_num, dr_dp):
            for nn, cc in zip(n, c):
                print(nn)
                print(cc)
                assert abs(nn - cc) < 1e-7

        ds_num = []
        for i in range(len(parameters)):

            def f(x):
                p = copy(parameters)
                p[i] = x
                return compute_sigma(p)

            fd = first_derivative(f, parameters[i], step)
            print(fd)
            ds_num.append(fd.reshape(3, 3, 1))
        ds_num = np.concatenate(ds_num, axis=2)

        for i in range(len(parameters)):
            for n, c in zip(ds_num[:, :, i], dS_dp[:, :, i]):
                for nn, cc in zip(n.flatten(), c.flatten()):
                    print(nn)
                    print(cc)
                    assert abs(nn - cc) < 1e-5

        dl_num = []
        for i in range(len(parameters)):

            def f(x):
                p = copy(parameters)
                p[i] = x
                return compute_sigma_lambda(p) ** 2

            dl_num.append(first_derivative(f, parameters[i], step))

        for n, c in zip(dl_num, dL_dp):
            assert abs(n - c) < 1e-7

    S1 = Simple1MosaicityParameterisation()
    S6 = Simple6MosaicityParameterisation()
    W = WavelengthSpreadParameterisation()

    check(S1, None, fix_wavelength_spread=True)
    check(S1, W, fix_mosaic_spread=True)
    check(S1, W, fix_wavelength_spread=True)
    check(S1, W, fix_unit_cell=True)
    check(S1, W, fix_orientation=True)

    check(S6, None, fix_wavelength_spread=True)
    check(S6, W, fix_mosaic_spread=True)
    check(S6, W, fix_wavelength_spread=True)
    check(S6, W, fix_unit_cell=True)
    check(S6, W, fix_orientation=True)
