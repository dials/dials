"""
Tests for the error model.
"""
from __future__ import absolute_import, division, print_function
from math import sqrt
import pytest
from libtbx.test_utils import approx_equal
from dials.algorithms.scaling.error_model.error_model import get_error_model
from dials.algorithms.scaling.error_model.error_model_target import ErrorModelTarget
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.scaling_refiner import error_model_refinery
from dials.array_family import flex
from cctbx.sgtbx import space_group


@pytest.fixture()
def large_reflection_table():
    """Create a reflection table."""
    return generate_refl_1()


@pytest.fixture(scope="module")
def test_sg():
    """Create a space group object."""
    return space_group("P 1")


def generate_refl_1():
    """Generate a test reflection table. Note tha the variance values are chosen
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


def data_for_error_model_test(background_variance=1, multiplicity=100, b=0.05):
    """Model a set of poisson-distributed observations on a constant-variance
    background."""

    ## First create a miller array of observations (in asu)
    from cctbx import miller
    from cctbx import crystal

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
    # sigmas need to be inflated by the error model with the given b.
    import scitbx
    from scitbx.random import variate, poisson_distribution

    scitbx.random.set_random_seed(0)
    intensities = flex.int()
    variances = flex.int()
    miller_index = flex.miller_index()
    for i, idx in zip(mean_intensities, ms.indices()):
        g = variate(poisson_distribution(mean=i))
        for _ in range(multiplicity):
            I = next(g)
            alpha = (1.0 + (b ** 2 * I)) ** 0.5
            intensities.append(int((alpha * I) + ((1.0 - alpha) * i)))
            variances.append(I + background_variance)
            miller_index.append(idx)

    reflections = flex.reflection_table()
    reflections["intensity"] = intensities.as_double()
    reflections["variance"] = variances.as_double()
    reflections["miller_index"] = miller_index
    reflections["inverse_scale_factor"] = flex.double(intensities.size(), 1.0)
    reflections["id"] = flex.int(intensities.size(), 1)

    return reflections


test_input = [
    [1, 100, [0.01, 2e-4], 0.0],  # no background var, no systematic error
    [5, 100, [0.01, 2e-4], 0.0],  # low background var
    [50, 100, [0.01, 2e-4], 0.0],  # high background var
    [1, 100, [0.01, 0.005], 0.05],  # no background var, signficant systematic error
    [5, 100, [0.01, 0.005], 0.05],  # low background var
    [50, 100, [0.01, 0.005], 0.05],  # high background var
    [1, 100, [0.01, 0.0025], 0.025],  # no background var, lower systematic error
    [5, 100, [0.01, 0.0025], 0.025],  # low background var
    [50, 100, [0.01, 0.0025], 0.025],  # high background var
    [50, 20, [0.02, 0.01], 0.05],  # low mult, high background, high system. error
    [50, 20, [0.02, 0.005], 0.025],  # low mult, high background, low systematic error
]
# Note, it appears that the b-parameter is usually slightly underestimated
# compared to the simulated data.


@pytest.mark.parametrize(
    "background_variance, multiplicity, abs_tolerances, model_b", test_input
)
def test_error_model_on_simulated_data(
    background_variance, multiplicity, abs_tolerances, model_b
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

    The purpose of the test is to show that the error model is finding the
    correct solution for a variety of background levels, with varying levels of
    systematic error (model b parameter, which is the one that has the highest
    effect on the errors output by scaling). Included are some more realistic
    cases, with lower multiplicity of measurement and high background variances.

    These tests are also designed to help validate the cutoff choices in the
    error model code:
    - the need for a min_Ih to ensure that poisson approx normal distribution
    - the need for an avg_I_over_var cutoff to remove background effects
    - cutting off extreme 'outlier' deviations as not to mislead the model."""

    data = data_for_error_model_test(
        int(background_variance), int(multiplicity), b=model_b
    )

    Ih_table = IhTable([data], space_group("P 2ac 2ab"))
    em = get_error_model("basic")
    block = Ih_table.blocked_data_list[0]
    em.min_reflections_required = 250
    error_model = em(block, n_bins=10)
    assert error_model.summation_matrix.n_rows > 400
    refinery = error_model_refinery(
        engine="SimpleLBFGS", target=ErrorModelTarget(error_model), max_iterations=100
    )
    refinery.run()
    error_model = refinery.return_error_model()
    assert error_model.refined_parameters[0] == pytest.approx(
        1.00, abs=abs_tolerances[0]
    )
    assert abs(error_model.refined_parameters[1]) == pytest.approx(
        model_b, abs=abs_tolerances[1]
    )


def test_errormodel(large_reflection_table, test_sg):
    """Test the initialisation and methods of the error model."""

    # first test get_error_model helper function.
    with pytest.raises(ValueError):
        em = get_error_model("bad")
    em = get_error_model("basic")
    em.min_reflections_required = 1
    Ih_table = IhTable([large_reflection_table], test_sg, nblocks=1)
    block = Ih_table.blocked_data_list[0]
    error_model = em(block, n_bins=2, min_Ih=1.0)
    assert error_model.summation_matrix[0, 1] == 1
    assert error_model.summation_matrix[1, 1] == 1
    assert error_model.summation_matrix[2, 0] == 1
    assert error_model.summation_matrix[3, 0] == 1
    assert error_model.summation_matrix[4, 0] == 1
    assert error_model.summation_matrix.non_zeroes == 5
    assert list(error_model.bin_counts) == [3, 2]

    # Test calc sigmaprime
    x0 = 1.0
    x1 = 0.1
    error_model.sigmaprime = error_model.calc_sigmaprime([x0, x1])
    cal_sigpr = list(
        x0
        * ((block.variances + ((x1 * block.intensities) ** 2)) ** 0.5)
        / block.inverse_scale_factors
    )
    assert list(error_model.sigmaprime) == cal_sigpr[4:7] + cal_sigpr[-2:]

    # Test calc delta_hl
    error_model.sigmaprime = error_model.calc_sigmaprime([1.0, 0.0])  # Reset
    # Calculate example for three elements, with intensities 1, 5 and 10 and
    # variances 1, 5 and 10 using he formula
    # delta_hl = sqrt(n_h - 1 / n_h) * (Ihl/ghl - Ih) / sigmaprime
    error_model.delta_hl = error_model.calc_deltahl()
    expected_deltas = [
        (-3.0 / 2.0) * sqrt(2.0 / 3.0),
        (5.0 / 2.0) * sqrt(2.0 / 15.0),
        5.0 * sqrt(2.0 / 30.0),
        -0.117647058824,
        0.124783549621,
    ]
    assert approx_equal(list(error_model.delta_hl), expected_deltas)


def test_error_model_target(large_reflection_table, test_sg):
    """Test the error model target."""
    Ih_table = IhTable([large_reflection_table], test_sg, nblocks=1)
    block = Ih_table.blocked_data_list[0]
    em = get_error_model("basic")
    em.min_reflections_required = 1
    error_model = em(block, n_bins=2, min_Ih=2.0)
    error_model.update_for_minimisation([1.0, 0.05])
    target = ErrorModelTarget(error_model, starting_values=[1.0, 0.05])
    # Test residual calculation
    residuals = target.calculate_residuals()
    assert residuals == (flex.double(2, 1.0) - error_model.bin_variances) ** 2

    # Test gradient calculation against finite differences.
    gradients = target.calculate_gradients()
    gradient_fd = calculate_gradient_fd(target)
    assert approx_equal(list(gradients), list(gradient_fd))

    # Test the method calls
    r, g = target.compute_functional_gradients()
    assert r == residuals
    assert list(gradients) == list(g)
    r, g, c = target.compute_functional_gradients_and_curvatures()
    assert r == residuals
    assert list(gradients) == list(g)
    assert c is None


def calculate_gradient_fd(target):
    """Calculate gradient array with finite difference approach."""
    delta = 1.0e-6
    gradients = flex.double([0.0] * len(target.x))
    # iterate over parameters, varying one at a time and calculating the gradient
    for i in range(len(target.x)):
        target.x[i] -= 0.5 * delta
        target.predict()
        R_low = target.calculate_residuals()
        target.x[i] += delta
        target.predict()
        R_upper = target.calculate_residuals()
        target.x[i] -= 0.5 * delta
        target.predict()
        gradients[i] = (flex.sum(R_upper) - flex.sum(R_low)) / delta
    return gradients
