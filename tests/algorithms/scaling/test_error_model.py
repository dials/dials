"""
Tests for the error model.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cctbx.sgtbx import space_group
from libtbx import phil

from dials.algorithms.scaling.error_model.engine import ErrorModelRefinery
from dials.algorithms.scaling.error_model.error_model import (
    BasicErrorModel,
    ErrorModelB_APM,
    calc_deltahl,
    calc_sigmaprime,
)
from dials.algorithms.scaling.error_model.error_model_target import ErrorModelTargetB
from dials.algorithms.scaling.Ih_table import IhTable
from dials.array_family import flex
from dials.util.options import ArgumentParser


@pytest.fixture()
def large_reflection_table():
    """Create a reflection table."""
    return generate_refl_1()


@pytest.fixture(scope="module")
def test_sg():
    """Create a space group object."""
    return space_group("P 1")


def generated_param():
    """Generate a param phil scope."""
    phil_scope = phil.parse(
        """
      include scope dials.command_line.scale.phil_scope
  """,
        process_includes=True,
    )
    parser = ArgumentParser(phil=phil_scope, check_format=False)
    parameters, _ = parser.parse_args(args=[], quick_parse=True, show_diff_phil=False)
    return parameters


def generate_refl_1():
    """Generate a test reflection table. Note that the variance values are chosen
    as the 'True' Ih_values, which would be found if unity weights were chosen
    in this example."""
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    )
    reflections["inverse_scale_factor"] = flex.double([1.0] * 5 + [2.0] * 5)
    reflections["variance"] = flex.double(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    )
    reflections["miller_index"] = flex.miller_index(
        [(1, 0, 0), (0, 0, 1), (-3, 0, 0), (0, 2, 0), (1, 0, 0)]
        + [(0, 0, -20), (0, 0, 2), (10, 0, 0), (10, 0, 0), (1, 0, 0)]
    )
    reflections.set_flags(flex.bool(10, True), reflections.flags.integrated)
    return reflections


def data_for_error_model_test(background_variance=1, multiplicity=100, b=0.05, a=1.0):
    """Model a set of poisson-distributed observations on a constant-variance
    background."""

    ## First create a miller array of observations (in asu)
    from cctbx import crystal, miller

    ms = miller.build_set(
        crystal_symmetry=crystal.symmetry(
            space_group_symbol="P212121", unit_cell=(12, 12, 25, 90, 90, 90)
        ),
        anomalous_flag=False,
        d_min=1.0,
    )
    assert ms.size() == 2150
    mean_intensities = 5.0 * (ms.d_spacings().data() ** 4)
    # ^ get a good range of intensities, with high intensity at low
    # miller index, mean = 285.2, median = 13.4

    # when applying b, use fact that I' - Imean = alpha(I - Imean), will
    # give the same distribution as sigma' = alpha sigma,
    # where alpha = (1 + (b^2 I)) ^ 0.5. i.e. this is the way to increase the
    # deviations of I-Imean and keep the same 'poisson' sigmas, such that the
    # sigmas need to be inflated by the error model with the given a, b.
    import scitbx
    from scitbx.random import poisson_distribution, variate

    # Note, if a and b both set, probably not quite right, but okay for small
    # a and b for the purpose of a test

    scitbx.random.set_random_seed(0)
    intensities = flex.int()
    variances = flex.double()
    miller_index = flex.miller_index()
    for i, idx in zip(mean_intensities, ms.indices()):
        g = variate(poisson_distribution(mean=i))
        for _ in range(multiplicity):
            intensity = next(g)
            if b > 0.0:
                alpha = (1.0 + (b**2 * intensity)) ** 0.5
                intensities.append(int((alpha * intensity) + ((1.0 - alpha) * i)))
            else:
                intensities.append(intensity)
            variances.append((intensity + background_variance) / (a**2))
            miller_index.append(idx)

    reflections = flex.reflection_table()
    reflections["intensity"] = intensities.as_double()
    reflections["variance"] = variances.as_double()
    reflections["miller_index"] = miller_index
    reflections["inverse_scale_factor"] = flex.double(intensities.size(), 1.0)
    reflections["id"] = flex.int(intensities.size(), 1)

    return reflections


no_bg_test_input = [
    # now background variance, high multiplicity
    [1, 100, [0.05, 5e-3], 1.0, 0.00],  # no background var, no systematic error
    [1, 100, [0.05, 5e-3], 1.0, 0.02],  # no background var, some systematic error
    [
        1,
        100,
        [0.05, 5e-3],
        1.0,
        0.04,
    ],  # no background var, significant systematic error
    [1, 100, [0.05, 5e-3], 1.4, 0.00],
]


@pytest.mark.parametrize(
    "background_variance, multiplicity, abs_tolerances, model_a, model_b",
    no_bg_test_input,
)
def test_error_model_on_simulated_data(
    background_variance, multiplicity, abs_tolerances, model_a, model_b
):
    """Test the refinement of the error model using simulated data.

    The simulated data consists of 2150 unique reflections, with I = 5.0 * d^4,
    giving an intensity range from 122070.3 to 5.0, with a mean of 285.21 and a
    median of 13.76.
    Each reflection is sampled 'multiplicity' times from a poisson distribution
    to give the intensities, which is scaled away from the mean I using the
    model_bm factor (scaled such that when the model is refined, the correct b
    should be returned).
    The test uses three different representative levels of background variance
    (1=None, 5=low, 50=high), which is added to the variance from poisson
    counting statistics. A tolerance of 10 % is used for the model b parameter
    for good cases, which is increased to 20 % for tough cases.
    """

    data = data_for_error_model_test(
        int(background_variance), int(multiplicity), b=model_b, a=model_a
    )

    Ih_table = IhTable([data], space_group("P 2ac 2ab"))

    block = Ih_table.blocked_data_list[0]
    BasicErrorModel.min_reflections_required = 250

    error_model = BasicErrorModel()
    error_model.configure_for_refinement(block)
    assert error_model.binner.summation_matrix.n_rows > 400
    refinery = ErrorModelRefinery(error_model, parameters_to_refine=["a", "b"])
    refinery.run()
    assert refinery.model.parameters[0] == pytest.approx(model_a, abs=abs_tolerances[0])
    assert abs(refinery.model.parameters[1]) == pytest.approx(
        model_b, abs=abs_tolerances[1]
    )


def test_errormodel(large_reflection_table, test_sg):
    """Test the initialisation and methods of the error model."""

    Ih_table = IhTable([large_reflection_table], test_sg, nblocks=1)
    block = Ih_table.blocked_data_list[0]
    params = generated_param()
    params.weighting.error_model.basic.n_bins = 2
    params.weighting.error_model.basic.min_Ih = 1.0
    em = BasicErrorModel
    em.min_reflections_required = 1
    error_model = em(basic_params=params.weighting.error_model.basic)
    error_model.min_reflections_required = 1
    error_model.configure_for_refinement(block)
    assert error_model.binner.summation_matrix[0, 1] == 1
    assert error_model.binner.summation_matrix[1, 1] == 1
    assert error_model.binner.summation_matrix[2, 0] == 1
    assert error_model.binner.summation_matrix[3, 0] == 1
    assert error_model.binner.summation_matrix[4, 0] == 1
    assert error_model.binner.summation_matrix.non_zeroes == 5
    assert list(error_model.binner.binning_info["refl_per_bin"]) == [3, 2]

    # Test calc sigmaprime
    x0 = 1.0
    x1 = 0.1
    sigmaprime = calc_sigmaprime([x0, x1], error_model.filtered_Ih_table)
    cal_sigpr = list(
        x0
        * np.sqrt(block.variances + np.square(x1 * block.intensities))
        / block.inverse_scale_factors
    )
    assert list(sigmaprime) == pytest.approx(cal_sigpr[4:7] + cal_sigpr[-2:])

    # Test calc delta_hl
    sigmaprime = calc_sigmaprime([1.0, 0.0], error_model.filtered_Ih_table)  # Reset
    # Calculate example for three elements, with intensities 1, 5 and 10 and
    # variances 1, 5 and 10 using he formula
    # delta_hl = math.sqrt(n_h - 1 / n_h) * (Ihl/ghl - Ih) / sigmaprime
    delta_hl = calc_deltahl(
        error_model.filtered_Ih_table,
        error_model.filtered_Ih_table.calc_nh(),
        sigmaprime,
    )
    expected_deltas = [
        (-3.0 / 2.0) * math.sqrt(2.0 / 3.0),
        (5.0 / 2.0) * math.sqrt(2.0 / 15.0),
        5.0 * math.sqrt(2.0 / 30.0),
        -0.117647058824,
        0.124783549621,
    ]
    assert list(delta_hl) == pytest.approx(expected_deltas)


def test_error_model_target(large_reflection_table, test_sg):
    """Test the error model target."""
    Ih_table = IhTable([large_reflection_table], test_sg, nblocks=1)
    block = Ih_table.blocked_data_list[0]
    em = BasicErrorModel
    em.min_reflections_required = 1
    params = generated_param()
    params.weighting.error_model.basic.n_bins = 2
    params.weighting.error_model.basic.min_Ih = 1.0
    error_model = em(basic_params=params.weighting.error_model.basic)
    error_model.configure_for_refinement(block)
    error_model.parameters = [1.0, 0.05]
    parameterisation = ErrorModelB_APM(error_model)
    target = ErrorModelTargetB(error_model)
    target.predict(parameterisation)

    # Test gradient calculation against finite differences.
    gradients = target.calculate_gradients(parameterisation)
    gradient_fd = calculate_gradient_fd(target, parameterisation)
    assert list(gradients) == pytest.approx(list(gradient_fd))

    # Test the method calls
    _, g = target.compute_functional_gradients(parameterisation)
    assert list(gradients) == pytest.approx(list(g))


def calculate_gradient_fd(target, parameterisation):
    """Calculate gradient array with finite difference approach."""
    delta = 1.0e-6
    parameterisation.set_param_vals([parameterisation.x[0] - (0.5 * delta)])
    target.predict(parameterisation)
    R_low = target.calculate_residuals(parameterisation)
    parameterisation.set_param_vals([parameterisation.x[0] + delta])
    target.predict(parameterisation)
    R_upper = target.calculate_residuals(parameterisation)
    parameterisation.set_param_vals([parameterisation.x[0] - (0.5 * delta)])
    target.predict(parameterisation)
    gradients = [(flex.sum(R_upper) - flex.sum(R_low)) / delta]
    return gradients
