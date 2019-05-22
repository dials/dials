from __future__ import division
from __future__ import print_function
from collections import namedtuple
from os.path import join
from math import exp
from random import uniform, randint
import pytest
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix
from dials.array_family import flex
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalUnitCellParameterisation,
    CrystalOrientationParameterisation,
)
from dials.algorithms.profile_model.potato.parameterisation import (
    Simple1MosaicityParameterisation,
    Simple6MosaicityParameterisation,
    # Angular2MosaicityParameterisation,
    # Angular4MosaicityParameterisation,
    # WavelengthSpreadParameterisation,
    ModelState,
    ReflectionModelState,
)
from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
    Simple1ProfileModel,
    Simple6ProfileModel,
)
from dials.algorithms.profile_model.potato.refiner import (
    ConditionalDistribution,
    rotate_vec3_double,
    rotate_mat3_double,
    ReflectionLikelihood,
    RefinerData,
    Refiner,
)


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
    mobs = matrix.col((uniform(0, 1e-3), uniform(0, 1e-3)))
    sp = s2.normalize() * s0.length()

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

    return params, s0, sp, h, ctot, mobs, Sobs


@pytest.fixture
def testdata(dials_regression):

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

    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)
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
                * exp(-0.5 * (j - 5) ** 2 / 1 ** 2)
                * exp(-0.5 * (i - 5) ** 2 / 1 ** 2)
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
        s0 = testdata.s0
        sp = testdata.sp
        h = testdata.h
        # ctot = testdata.ctot
        # mobs = testdata.mobs
        # Sobs = testdata.Sobs

        U_params = models[1].get_param_vals()
        B_params = models[2].get_param_vals()
        M_params = flex.double(models[0][: mosaicity_parameterisation.num_parameters()])
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
        state.set_U_params(U_params)
        state.set_B_params(B_params)
        state.set_M_params(M_params)
        state.set_L_params(L_params)

        model = ReflectionModelState(state, s0, h)

        def get_conditional(model):
            # Compute the change of basis
            R = compute_change_of_basis_operation(s0, sp)

            # The s2 vector
            r = model.get_r()
            s2 = s0 + r

            # Rotate the mean vector
            mu = R * s2

            # Rotate the covariance matrix
            S = R * model.get_sigma() * R.transpose()

            # Rotate the first derivative matrices
            dS = rotate_mat3_double(R, model.get_dS_dp())

            # Rotate the first derivative of s2
            dmu = rotate_vec3_double(R, model.get_dr_dp())

            # Construct the conditional distribution
            conditional = ConditionalDistribution(s0, mu, dmu, S, dS)
            return conditional

        conditional = get_conditional(model)

        step = 1e-6

        dm_dp = conditional.first_derivatives_of_mean()
        dS_dp = conditional.first_derivatives_of_sigma()

        parameters = state.get_active_parameters()

        def compute_sigma(parameters):
            state.set_active_parameters(parameters)
            model = ReflectionModelState(state, s0, h)
            conditional = get_conditional(model)
            return conditional.sigma()

        def compute_mean(parameters):
            state.set_active_parameters(parameters)
            model = ReflectionModelState(state, s0, h)
            conditional = get_conditional(model)
            return conditional.mean()

        dm_num = []
        for i in range(len(parameters)):

            def f(x):
                p = [pp for pp in parameters]
                p[i] = x
                return compute_mean(p)

            dm_num.append(first_derivative(f, parameters[i], step))

        for n, c in zip(dm_num, dm_dp):
            assert all(abs(nn - cc) < 1e-7 for nn, cc in zip(n, c))

        ds_num = []
        for i in range(len(parameters)):

            def f(x):
                p = [pp for pp in parameters]
                p[i] = x
                return compute_sigma(p)

            ds_num.append(first_derivative(f, parameters[i], step))

        for n, c in zip(ds_num, dS_dp):
            assert all(abs(nn - cc) < 1e-7 for nn, cc in zip(n, c))

    S1 = Simple1MosaicityParameterisation()
    S6 = Simple6MosaicityParameterisation()

    check(S1, None, fix_wavelength_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S1, None, fix_wavelength_spread=True, fix_orientation=True)

    check(S6, None, fix_wavelength_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S6, None, fix_wavelength_spread=True, fix_orientation=True)


def test_rotate_vec3_double():

    vectors = flex.vec3_double([matrix.col((1, 1, 1)).normalize()])

    R = compute_change_of_basis_operation(matrix.col((0, 0, 1)), matrix.col(vectors[0]))

    rotated = rotate_vec3_double(R, vectors)

    assert rotated[0] == pytest.approx((0, 0, 1))


def test_rotate_mat3_double():

    A = matrix.diag((1, 1, 1))
    R = compute_change_of_basis_operation(matrix.col((0, 0, 1)), matrix.col((1, 1, 1)))
    A = R.transpose() * A * R
    matrices = flex.mat3_double([A])

    R = R.transpose()

    rotated = rotate_mat3_double(R, matrices)

    assert rotated[0] == pytest.approx((1, 0, 0, 0, 1, 0, 0, 0, 1))


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
        s0 = testdata.s0
        sp = testdata.sp
        h = testdata.h
        ctot = testdata.ctot
        mobs = testdata.mobs
        Sobs = testdata.Sobs

        U_params = models[1].get_param_vals()
        B_params = models[2].get_param_vals()
        M_params = flex.double(models[0][: mosaicity_parameterisation.num_parameters()])
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
        state.set_U_params(U_params)
        state.set_B_params(B_params)
        state.set_M_params(M_params)
        state.set_L_params(L_params)

        def get_reflection_likelihood(state):
            return ReflectionLikelihood(state, s0, sp, h, ctot, mobs, Sobs)

        likelihood = get_reflection_likelihood(state)

        step = 1e-6

        dL_dp = likelihood.first_derivatives()

        parameters = state.get_active_parameters()

        assert len(dL_dp) == len(parameters)

        def compute_likelihood(parameters):
            state.set_active_parameters(parameters)
            likelihood = get_reflection_likelihood(state)
            return likelihood.log_likelihood()

        dL_num = []
        for i in range(len(parameters)):

            def f(x):
                p = [pp for pp in parameters]
                p[i] = x
                return compute_likelihood(p)

            dL_num.append(first_derivative(f, parameters[i], step))

        assert len(dL_num) == len(parameters)
        for n, c in zip(dL_num, dL_dp):
            assert n == pytest.approx(c)

    S1 = Simple1MosaicityParameterisation()
    S6 = Simple6MosaicityParameterisation()

    check(S1, None, fix_wavelength_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S1, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S1, None, fix_wavelength_spread=True, fix_orientation=True)

    check(S6, None, fix_wavelength_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_mosaic_spread=True)
    check(S6, None, fix_wavelength_spread=True, fix_unit_cell=True)
    check(S6, None, fix_wavelength_spread=True, fix_orientation=True)


def test_Refiner(testdata, refinerdata_testdata):

    experiment = testdata.experiment
    data = refinerdata_testdata

    sigma_d = 0.02 ** 2
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

        # print(state.get_unit_cell().parameters())
        print(state.get_M())

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
                * exp(-0.5 * (j - 5) ** 2 / 1 ** 2)
                * exp(-0.5 * (i - 5) ** 2 / 1 ** 2)
            )
            shoebox_mask[0, j, i] = 5
    for sbox in reflections["shoebox"]:
        sbox.data = shoebox_data
        sbox.mask = shoebox_mask

    data = RefinerData.from_reflections(experiment, reflections)

    assert tuple(data.s0) == pytest.approx(experiment.beam.get_s0())
    assert data.h_list == reflections["miller_index"]
    for a, b in zip(data.sp_list, reflections["sp"]):
        assert a == pytest.approx(b)
    assert data.ctot_list == sum(shoebox_data)
    mobs1, mobs2 = data.mobs_list.parts()
    mobs1 = flex.abs(mobs1)
    mobs2 = flex.abs(mobs2)
    assert max(mobs1) < 1e-6
    assert max(mobs2) < 1e-6
