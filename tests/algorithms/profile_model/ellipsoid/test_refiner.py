from __future__ import annotations

import copy
from collections import namedtuple
from math import exp
from random import randint, uniform

import numpy as np
import pytest
from numpy.linalg import norm

from scitbx import matrix

from dials.algorithms.profile_model.ellipsoid.model import (
    Simple1ProfileModel,
    Simple6ProfileModel,
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.ellipsoid.parameterisation import (
    ModelState,
    Simple1Angular1MosaicityParameterisation,
    Simple1Angular3MosaicityParameterisation,
    Simple1MosaicityParameterisation,
    Simple6Angular1MosaicityParameterisation,
    Simple6Angular3MosaicityParameterisation,
    Simple6MosaicityParameterisation,
)
from dials.algorithms.profile_model.ellipsoid.refiner import (
    Refiner,
    RefinerData,
    ReflectionLikelihood,
    rotate_mat3_double,
    rotate_vec3_double,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)
from dials.array_family import flex


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
    mobs = np.array([uniform(0, 1e-3), uniform(0, 1e-3)])
    sp = s2.normalize() * s0.length()

    b1, b2, b3, b4, b5, b6, b7, b8, b9 = (
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(0.0, 1e-3),
        uniform(0.0, 1e-3),
    )

    S_param = (b1, b2, b3, b4, b5, b6, b7, b8, b9)
    L_param = (uniform(1e-3, 2e-3),)
    ctot = randint(100, 1000)

    T = np.array(
        (uniform(1e-3, 2e-3), 0, uniform(1e-6, 2e-6), uniform(1e-3, 2e-3))
    ).reshape(2, 2)
    Sobs = np.matmul(T, T.T)

    params = [S_param, U_param, B_param, L_param]

    return params, s0, sp, h, ctot, mobs, Sobs


@pytest.fixture
def testdata(test_experiment):

    TestData = namedtuple(
        "TestData",
        [
            "experiment",
            "reflections",
            "models",
            "s0",
            "sp",
            "h",
            "ctot",
            "mobs",
            "Sobs",
        ],
    )

    experiments = [test_experiment]
    reflections = flex.reflection_table.from_predictions_multi(experiments)

    models, s0, sp, h, ctot, mobs, Sobs = generate_data(experiments, reflections)

    return TestData(
        experiment=experiments[0],
        reflections=reflections,
        models=models,
        s0=s0,
        sp=sp,
        h=h,
        ctot=ctot,
        mobs=mobs,
        Sobs=Sobs,
    )


@pytest.fixture
def refinerdata_testdata(testdata):

    experiment = testdata.experiment
    reflections = testdata.reflections

    panel = experiment.detector[0]
    s0_length = matrix.col(experiment.beam.get_s0()).length()
    reflections["bbox"] = flex.int6(len(reflections))
    reflections["xyzobs.px.value"] = flex.vec3_double(len(reflections))
    reflections["s2"] = reflections["s1"].each_normalize() * s0_length
    reflections["sp"] = flex.vec3_double(len(reflections))
    for i, (x, y, z) in enumerate(reflections["xyzcal.px"]):
        x0 = int(x) - 5
        x1 = int(x) + 5 + 1
        y0 = int(y) - 5
        y1 = int(y) + 5 + 1
        z0 = int(z)
        z1 = z0 + 1
        reflections["bbox"][i] = x0, x1, y0, y1, z0, z1
        reflections["xyzobs.px.value"][i] = (int(x) + 0.5, int(y) + 0.5, int(z) + 0.5)
        reflections["sp"][i] = (
            matrix.col(
                panel.get_pixel_lab_coord(reflections["xyzobs.px.value"][i][0:2])
            ).normalize()
            * s0_length
        )

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"], reflections["bbox"], allocate=True
    )

    shoebox_data = flex.float(flex.grid(1, 11, 11))
    shoebox_mask = flex.int(flex.grid(1, 11, 11))
    for j in range(11):
        for i in range(11):
            shoebox_data[0, j, i] = (
                100
                * exp(-0.5 * (j - 5) ** 2 / 1**2)
                * exp(-0.5 * (i - 5) ** 2 / 1**2)
            )
            shoebox_mask[0, j, i] = 5
    for sbox in reflections["shoebox"]:
        sbox.data = shoebox_data
        sbox.mask = shoebox_mask

    return RefinerData.from_reflections(experiment, reflections)


def test_ConditionalDistribution(testdata):
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
        h = testdata.h
        s0 = np.array(testdata.s0).reshape(3, 1)
        sp = np.array(testdata.sp).reshape(3, 1)
        ctot = testdata.ctot
        mobs = testdata.mobs
        Sobs = testdata.Sobs

        U_params = models[1].get_param_vals()
        B_params = models[2].get_param_vals()
        M_params = models[0][: mosaicity_parameterisation.num_parameters()]
        L_params = models[3]

        state = ModelState(
            experiment,
            mosaicity_parameterisation,
            wavelength_parameterisation,
            fix_mosaic_spread=fix_mosaic_spread,
            fix_wavelength_spread=fix_wavelength_spread,
            fix_unit_cell=fix_unit_cell,
            fix_orientation=fix_orientation,
        )
        state.U_params = np.array(U_params, dtype=np.float64)
        state.B_params = np.array(B_params, dtype=np.float64)
        state.M_params = np.array(M_params, dtype=np.float64)
        state.L_params = L_params

        def get_conditional(state):
            conditional = ReflectionLikelihood(
                state,
                s0,
                sp,
                h,
                ctot,
                np.array(mobs),
                np.array(Sobs).reshape(2, 2),
            ).conditional
            return conditional

        conditional = get_conditional(state)

        step = 1e-6

        dm_dp = copy.copy(conditional.first_derivatives_of_mean())
        dS_dp = copy.copy(conditional.first_derivatives_of_sigma())

        parameters = state.active_parameters

        def compute_sigma(parameters):
            state.active_parameters = parameters
            return get_conditional(state).sigma()

        def compute_mean(parameters):
            state.active_parameters = parameters
            return get_conditional(state).mean()

        dm_num = []
        for i in range(len(parameters)):

            def f(x):
                p = copy.copy(parameters)
                p[i] = x
                return compute_mean(p)

            dm_num.append(first_derivative(f, parameters[i], step))

        for n, c in zip(dm_num, dm_dp):
            for nn, cc in zip(n, c):
                if not abs(nn - cc) < 1e-7:
                    print(dm_num)
                    print(dm_dp)
                    print(nn)
                    print(cc)
                    assert 0
            # assert all(abs(nn - cc) < 1e-7 for nn, cc in zip(n, c))

        ds_num = []
        for i in range(len(parameters)):

            def f(x):
                p = copy.copy(parameters)
                p[i] = x
                return compute_sigma(p)

            ds_num.append(first_derivative(f, parameters[i], step))

        for n, c in zip(ds_num, dS_dp):
            for nn, cc in zip(n.flatten(), c.flatten()):
                if not abs(nn - cc) < 1e-7:
                    print(nn)
                    print(cc)
                    assert 0
            # assert all(abs(nn - cc) < 1e-7 for nn, cc in zip(n, c))

    S1 = Simple1MosaicityParameterisation(np.array([0.01]))
    S6 = Simple6MosaicityParameterisation(
        np.array([0.01, 0.005, 0.02, 0.015, 0.03, 0.025])
    )
    S1A1 = Simple1Angular1MosaicityParameterisation(np.array([0.01, 0.002]))
    S6A1 = Simple6Angular1MosaicityParameterisation(
        np.array([0.01, 0.005, 0.02, 0.015, 0.03, 0.025, 0.002])
    )
    S1A3 = Simple1Angular3MosaicityParameterisation(
        np.array([0.01, 0.002, 0.001, 0.002])
    )
    S6A3 = Simple6Angular3MosaicityParameterisation(
        np.array([0.01, 0.005, 0.02, 0.015, 0.03, 0.025, 0.002, 0.001, 0.003])
    )

    check(S1, None, fix_wavelength_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S1, None, fix_wavelength_spread=True, fix_orientation=True)

    check(S6, None, fix_wavelength_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S6, None, fix_wavelength_spread=True, fix_orientation=True)

    # Note, we cannot use finite difference tests to test the derivatives when
    # the unit cell and orientation are not fixed for the angular model. This is because
    # the derivative calculation is under the assumption that the Q matrix is invariant,
    # however the Q matrix depends on r so there is a small additional term missing (dQ/dp).
    # So just test the mosaic model minimisation at fixed UB.
    check(
        S1A1, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )
    check(
        S6A1, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )
    check(
        S1A3, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )
    check(
        S6A3, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )


def test_rotate_vec3_double():

    vectors = np.array([[1, 1, 1]], dtype=np.float64).reshape(3, 1) / norm(
        np.array([1, 1, 1], dtype=np.float64)
    )

    v1 = np.array([0, 0, 1.0], dtype=np.float64).reshape(3, 1)
    R = compute_change_of_basis_operation(v1, vectors)

    rotated = rotate_vec3_double(R, vectors)

    assert rotated == pytest.approx((0, 0, 1))


def test_rotate_mat3_double():

    A = np.eye(3, dtype=np.float64)
    v1 = np.array([0, 0, 1.0], dtype=np.float64).reshape(3, 1)
    v2 = np.array([1, 1, 1.0], dtype=np.float64).reshape(3, 1)
    R = compute_change_of_basis_operation(v1, v2)
    A = np.matmul(np.matmul(R.T, A), R).reshape(3, 3, 1)

    R = R.T

    rotated = rotate_mat3_double(R, A)

    assert rotated[:, :, 0].flatten() == pytest.approx(
        (1, 0, 0, 0, 1, 0, 0, 0, 1), abs=1e-12
    )


def test_ReflectionLikelihood(testdata):
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
        s0 = np.array(testdata.s0)
        sp = np.array(testdata.sp)
        h = testdata.h
        ctot = testdata.ctot
        mobs = testdata.mobs
        Sobs = testdata.Sobs

        U_params = models[1].get_param_vals()
        B_params = models[2].get_param_vals()
        M_params = np.array(models[0][: mosaicity_parameterisation.num_parameters()])
        L_params = flex.double(models[3])

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

        def get_reflection_likelihood(state):
            return ReflectionLikelihood(state, s0, sp, h, ctot, mobs, Sobs)

        likelihood = get_reflection_likelihood(state)

        step = 1e-6

        dL_dp = likelihood.first_derivatives()

        parameters = state.active_parameters

        assert len(dL_dp) == len(parameters)

        def compute_likelihood(parameters):
            state.active_parameters = parameters
            likelihood = get_reflection_likelihood(state)
            return likelihood.log_likelihood()

        dL_num = []
        for i in range(len(parameters)):

            def f(x):
                p = copy.copy(parameters)
                p[i] = x
                return compute_likelihood(p)

            dL_num.append(first_derivative(f, parameters[i], step))

        assert len(dL_num) == len(parameters)
        print(dL_num)
        print(list(dL_dp))
        for n, c in zip(dL_num, dL_dp):
            print(n, c)
            assert n == pytest.approx(c, rel=1e-6)

    S1 = Simple1MosaicityParameterisation()
    S6 = Simple6MosaicityParameterisation()
    S1A1 = Simple1Angular1MosaicityParameterisation(np.array([0.01, 0.002]))
    S6A1 = Simple6Angular1MosaicityParameterisation(
        np.array([0.01, 0.005, 0.02, 0.015, 0.03, 0.025, 0.002])
    )
    S1A3 = Simple1Angular3MosaicityParameterisation(
        np.array([0.01, 0.002, 0.001, 0.003])
    )
    S6A3 = Simple6Angular3MosaicityParameterisation(
        np.array([0.01, 0.005, 0.02, 0.015, 0.03, 0.025, 0.002, 0.001, 0.003])
    )

    check(S1, None, fix_wavelength_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S1, None, fix_wavelength_spread=True, fix_orientation=True)

    check(S6, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S6, None, fix_wavelength_spread=True, fix_orientation=True)
    check(S6, None, fix_wavelength_spread=True)

    # See note at L300 regarding finite difference tests at fixed UB.
    check(
        S1A1, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )
    check(
        S6A1, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )
    check(
        S1A3, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )
    check(
        S6A3, None, fix_wavelength_spread=True, fix_unit_cell=True, fix_orientation=True
    )


def test_Refiner(testdata, refinerdata_testdata):

    experiment = testdata.experiment
    data = refinerdata_testdata

    sigma_d = 0.02**2
    S1 = Simple1ProfileModel.from_sigma_d(sigma_d).parameterisation()
    S6 = Simple6ProfileModel.from_sigma_d(sigma_d).parameterisation()

    def check(
        parameterisation,
        fix_mosaic_spread=True,
        fix_orientation=True,
        fix_unit_cell=True,
        fix_wavelength_spread=True,
    ):

        state = ModelState(
            experiment,
            parameterisation,
            fix_mosaic_spread=fix_mosaic_spread,
            fix_orientation=fix_orientation,
            fix_unit_cell=fix_unit_cell,
            fix_wavelength_spread=fix_wavelength_spread,
        )

        refiner = Refiner(state, data)
        refiner.refine()

    check(S1, fix_mosaic_spread=False, fix_orientation=True, fix_unit_cell=True)
    check(S1, fix_mosaic_spread=True, fix_orientation=False, fix_unit_cell=True)
    check(S1, fix_mosaic_spread=True, fix_orientation=True, fix_unit_cell=False)
    check(S1, fix_mosaic_spread=False, fix_orientation=False, fix_unit_cell=True)
    check(S1, fix_mosaic_spread=True, fix_orientation=False, fix_unit_cell=False)
    check(S1, fix_mosaic_spread=False, fix_orientation=False, fix_unit_cell=False)

    check(S6, fix_mosaic_spread=False, fix_orientation=True, fix_unit_cell=True)
    check(S6, fix_mosaic_spread=True, fix_orientation=False, fix_unit_cell=True)
    check(S6, fix_mosaic_spread=True, fix_orientation=True, fix_unit_cell=False)
    check(S6, fix_mosaic_spread=False, fix_orientation=False, fix_unit_cell=True)
    check(S6, fix_mosaic_spread=True, fix_orientation=False, fix_unit_cell=False)
    check(S6, fix_mosaic_spread=False, fix_orientation=False, fix_unit_cell=False)


def test_RefinerData(testdata):

    experiment = testdata.experiment
    reflections = testdata.reflections

    panel = experiment.detector[0]
    s0_length = matrix.col(experiment.beam.get_s0()).length()
    reflections["bbox"] = flex.int6(len(reflections))
    reflections["xyzobs.px.value"] = flex.vec3_double(len(reflections))
    reflections["s2"] = reflections["s1"].each_normalize() * s0_length
    reflections["sp"] = flex.vec3_double(len(reflections))
    for i, (x, y, z) in enumerate(reflections["xyzcal.px"]):
        x0 = int(x) - 5
        x1 = int(x) + 5 + 1
        y0 = int(y) - 5
        y1 = int(y) + 5 + 1
        z0 = int(z)
        z1 = z0 + 1
        reflections["bbox"][i] = x0, x1, y0, y1, z0, z1
        reflections["xyzobs.px.value"][i] = (int(x) + 0.5, int(y) + 0.5, int(z) + 0.5)
        reflections["sp"][i] = (
            matrix.col(
                panel.get_pixel_lab_coord(reflections["xyzobs.px.value"][i][0:2])
            ).normalize()
            * s0_length
        )

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"], reflections["bbox"], allocate=True
    )

    shoebox_data = flex.float(flex.grid(1, 11, 11))
    shoebox_mask = flex.int(flex.grid(1, 11, 11))
    for j in range(11):
        for i in range(11):
            shoebox_data[0, j, i] = (
                100
                * exp(-0.5 * (j - 5) ** 2 / 1**2)
                * exp(-0.5 * (i - 5) ** 2 / 1**2)
            )
            shoebox_mask[0, j, i] = 5
    for sbox in reflections["shoebox"]:
        sbox.data = shoebox_data
        sbox.mask = shoebox_mask

    data = RefinerData.from_reflections(experiment, reflections)

    assert tuple(data.s0) == pytest.approx(experiment.beam.get_s0())
    assert data.h_list == reflections["miller_index"]
    for i, sp in enumerate(reflections["sp"]):
        assert data.sp_list[:, i] == pytest.approx(sp)
    assert data.ctot_list[0] == sum(shoebox_data)

    mobs1 = np.abs(data.mobs_list[0, :])
    mobs2 = np.abs(data.mobs_list[1, :])
    assert np.max(mobs1) < 1e-6
    assert np.max(mobs2) < 1e-6
