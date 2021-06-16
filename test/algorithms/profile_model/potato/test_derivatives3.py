from __future__ import division, print_function

from collections import namedtuple
from os.path import join
from random import randint, uniform

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix

from dials.algorithms.profile_model.potato.parameterisation import (
    ModelState,
    Simple6MosaicityParameterisation,
)
from dials.algorithms.profile_model.potato.refiner import ReflectionLikelihood
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

    b1, b2, b3, b4, b5, b6 = (
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
    )

    params = (b1, b2, b3, b4, b5, b6)

    S_param = Simple6MosaicityParameterisation(params)
    L_param = (uniform(1e-3, 2e-3),)
    ctot = randint(100, 1000)

    T = matrix.sqr((uniform(1e-3, 2e-3), 0, uniform(1e-6, 2e-6), uniform(1e-3, 2e-3)))
    Sobs = T * T.transpose()

    params = [S_param, U_param, B_param, L_param]

    return params, s0, h, ctot, mobs, Sobs


@pytest.fixture
def testdata(dials_regression):

    TestData = namedtuple(
        "TestData", ["experiment", "models", "s0", "h", "ctot", "mobs", "Sobs"]
    )

    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)
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


def test_first_derivatives(testdata):
    # this was called tmptest - complete?

    experiment = testdata.experiment
    models = testdata.models
    s0 = testdata.s0
    h = testdata.h
    ctot = testdata.ctot
    mobs = testdata.mobs
    Sobs = testdata.Sobs

    U_params = models[1].get_param_vals()
    B_params = models[2].get_param_vals()
    M_params = flex.double(models[0].parameters())

    L_params = flex.double(models[3])

    parameterisation = Simple6MosaicityParameterisation()
    state = ModelState(experiment, parameterisation)
    state.set_U_params(U_params)
    state.set_B_params(B_params)
    state.set_M_params(M_params)
    state.set_L_params(L_params)

    def compute_L(parameters):
        state.set_active_parameters(parameters)
        rd = ReflectionLikelihood(state, s0, mobs, h, ctot, matrix.col((0, 0)), Sobs)
        return rd.log_likelihood()

    step = 1e-5

    parameters = state.get_active_parameters()

    dL_num = []
    for i in range(len(parameters)):

        def f(x):
            p = [pp for pp in parameters]
            p[i] = x
            return compute_L(p)

        dL_num.append(first_derivative(f, parameters[i], step))

    state.set_active_parameters(parameters)
    rd = ReflectionLikelihood(state, s0, mobs, h, ctot, matrix.col((0, 0)), Sobs)

    dL_cal = list(rd.first_derivatives())

    assert all(abs(n - c) < 1e-3 for n, c in zip(dL_num, dL_cal))
