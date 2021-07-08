from math import sqrt
from os.path import join

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix

from dials.algorithms.profile_model.potato import chisq_quantile
from dials.algorithms.profile_model.potato.model import (
    Simple1ProfileModel,
    Simple6ProfileModel,
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.potato.parameterisation import (
    ModelState,
    Simple1MosaicityParameterisation,
    Simple6MosaicityParameterisation,
)
from dials.algorithms.spot_prediction import IndexGenerator
from dials.array_family import flex


@pytest.fixture
def simple1_profile_model():

    params = flex.double([4e-4])

    model = Simple1ProfileModel(params)
    return model


@pytest.fixture
def simple6_profile_model():
    params = flex.double([4e-4, 1e-4, 5e-4, 2e-4, 3e-4, 6e-4])
    model = Simple6ProfileModel(params)
    return model


@pytest.fixture
def simple1_model_state(dials_regression):

    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    state = ModelState(experiments[0], Simple1MosaicityParameterisation())

    return state


@pytest.fixture
def simple6_model_state(dials_regression):

    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    state = ModelState(experiments[0], Simple6MosaicityParameterisation())

    return state


def check_simple1_sigma(sigma, params):
    b1 = params[0]

    assert sigma[0] == pytest.approx(b1 ** 2)
    assert sigma[1] == 0
    assert sigma[2] == 0
    assert sigma[3] == 0
    assert sigma[4] == pytest.approx(b1 ** 2)
    assert sigma[5] == 0
    assert sigma[6] == 0
    assert sigma[7] == 0
    assert sigma[8] == pytest.approx(b1 ** 2)


def check_simple6_sigma(sigma, params):

    b1, b2, b3, b4, b5, b6 = params
    L = matrix.sqr((b1, 0, 0, b2, b3, 0, b4, b5, b6))
    M = L * L.transpose()
    for i in range(9):
        assert sigma[i] == pytest.approx(M[i])


def test_Simple1ProfileModel_sigma(simple1_profile_model):
    sigma = simple1_profile_model.sigma()
    check_simple1_sigma(sigma, [4e-4])


def test_Simple1ProfileModel_update_model(simple1_profile_model, simple1_model_state):
    params = [5e-4]
    simple1_model_state.set_M_params(params)
    simple1_profile_model.update_model(simple1_model_state)
    sigma = simple1_profile_model.sigma()
    check_simple1_sigma(sigma, params)


def test_Simple1ProfileModel_predict_reflections(
    dials_regression, simple1_profile_model
):

    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        2.0,
    )

    # Get an array of miller indices
    miller_indices = index_generator.to_array()
    reflections = simple1_profile_model.predict_reflections(
        experiments, miller_indices, probability=0.9973
    )

    # s2 = reflections["s2"]
    s0 = matrix.col(experiments[0].beam.get_s0())
    quantile = chisq_quantile(3, 0.9973)
    sigma_inv = matrix.sqr(simple1_profile_model.sigma()).inverse()

    for s2 in reflections["s2"]:
        s2_ = matrix.col(s2)
        x = s2_.normalize() * s0.length() - s2_
        d = (x.transpose() * sigma_inv * x)[0]
        assert d < quantile


def test_Simple1ProfileModel_compute_bbox(dials_regression, simple1_profile_model):
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        2.0,
    )

    # Get an array of miller indices
    miller_indices = index_generator.to_array()
    reflections = simple1_profile_model.predict_reflections(
        experiments, miller_indices, probability=0.9973
    )

    simple1_profile_model.compute_bbox(experiments, reflections)


def test_Simple1ProfileModel_compute_mask(dials_regression, simple1_profile_model):
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        2.0,
    )

    # Get an array of miller indices
    miller_indices = index_generator.to_array()
    reflections = simple1_profile_model.predict_reflections(
        experiments, miller_indices, probability=0.9973
    )

    simple1_profile_model.compute_bbox(experiments, reflections)

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"], reflections["bbox"], allocate=True
    )

    simple1_profile_model.compute_mask(experiments, reflections)


def test_Simple1ProfileModel_sigma_for_reflection(simple1_profile_model):
    sigma = simple1_profile_model.sigma()
    check_simple1_sigma(sigma, [4e-4])


def test_Simple1ProfileModel_parameterisation(simple1_profile_model):
    assert isinstance(
        simple1_profile_model.parameterisation(), Simple1MosaicityParameterisation
    )


def test_Simple6ProfileModel_sigma(simple6_profile_model):
    sigma = simple6_profile_model.sigma()
    check_simple6_sigma(sigma, [4e-4, 1e-4, 5e-4, 2e-4, 3e-4, 6e-4])


def test_Simple6ProfileModel_update_model(simple6_profile_model, simple6_model_state):
    params = [5e-4, 2e-4, 6e-4, 3e-4, 4e-4, 7e-4]
    simple6_model_state.set_M_params(params)
    simple6_profile_model.update_model(simple6_model_state)
    sigma = simple6_profile_model.sigma()
    check_simple6_sigma(sigma, params)


def test_Simple6ProfileModel_predict_reflections(
    dials_regression, simple6_profile_model
):
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        2.0,
    )

    # Get an array of miller indices
    miller_indices = index_generator.to_array()
    reflections = simple6_profile_model.predict_reflections(
        experiments, miller_indices, probability=0.9973
    )

    s2 = reflections["s2"]
    s0 = matrix.col(experiments[0].beam.get_s0())
    quantile = chisq_quantile(3, 0.9973)
    sigma_inv = matrix.sqr(simple6_profile_model.sigma()).inverse()

    for s2 in map(matrix.col, reflections["s2"]):
        x = s2.normalize() * s0.length() - s2
        d = (x.transpose() * sigma_inv * x)[0]
        assert d < quantile


def test_Simple6ProfileModel_compute_bbox(dials_regression, simple6_profile_model):
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        2.0,
    )

    # Get an array of miller indices
    miller_indices = index_generator.to_array()
    reflections = simple6_profile_model.predict_reflections(
        experiments, miller_indices, probability=0.9973
    )

    s2 = reflections["s2"]
    s0 = matrix.col(experiments[0].beam.get_s0())
    quantile = chisq_quantile(3, 0.9973)
    sigma_inv = matrix.sqr(simple6_profile_model.sigma()).inverse()

    for s2 in map(matrix.col, reflections["s2"]):
        x = s2.normalize() * s0.length() - s2
        d = (x.transpose() * sigma_inv * x)[0]
        assert d < quantile

    simple6_profile_model.compute_bbox(experiments, reflections)


def test_Simple6ProfileModel_compute_mask(dials_regression, simple6_profile_model):
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    # Create the index generator
    index_generator = IndexGenerator(
        experiments[0].crystal.get_unit_cell(),
        experiments[0].crystal.get_space_group().type(),
        2.0,
    )

    # Get an array of miller indices
    miller_indices = index_generator.to_array()
    reflections = simple6_profile_model.predict_reflections(
        experiments, miller_indices, probability=0.9973
    )

    s2 = reflections["s2"]
    s0 = matrix.col(experiments[0].beam.get_s0())
    quantile = chisq_quantile(3, 0.9973)
    sigma_inv = matrix.sqr(simple6_profile_model.sigma()).inverse()

    for s2 in map(matrix.col, reflections["s2"]):
        x = s2.normalize() * s0.length() - s2
        d = (x.transpose() * sigma_inv * x)[0]
        assert d < quantile

    simple6_profile_model.compute_bbox(experiments, reflections)

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"], reflections["bbox"], allocate=True
    )

    simple6_profile_model.compute_mask(experiments, reflections)


def test_Simple6ProfileModel_sigma_for_reflection(simple6_profile_model):
    sigma = simple6_profile_model.sigma()
    check_simple6_sigma(sigma, [4e-4, 1e-4, 5e-4, 2e-4, 3e-4, 6e-4])


def test_Simple6ProfileModel_parameterisation(simple6_profile_model):
    assert isinstance(
        simple6_profile_model.parameterisation(), Simple6MosaicityParameterisation
    )


def test_compute_change_of_basis_operation():

    r = matrix.col((0, 0.5, -(1 - sqrt(0.75))))
    s0 = matrix.col((0, 0, 1))
    s2 = s0 + r

    R = compute_change_of_basis_operation(s0, s2)

    z = R * s2

    assert z.length() == pytest.approx(1)
