from __future__ import division, print_function

from math import log, sqrt

import numpy.random
import pytest

from scitbx import matrix

from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.potato.parameterisation import (
    Simple6MosaicityParameterisation,
)
from dials.algorithms.profile_model.potato.refiner import MaximumLikelihoodTarget
from dials.algorithms.profile_model.potato.refiner import Refiner as ProfileRefiner
from dials.algorithms.profile_model.potato.refiner import (
    RefinerData as ProfileRefinerData,
)
from dials.algorithms.profile_model.potato.util.generate_simple import (
    generate_simple,
    generate_simple_binned,
)
from dials.algorithms.profile_model.potato.util.simplex import SimpleSimplex
from dials.array_family import flex


def log_likelihood(params, s0, s2_list, xbar_list, ctot_list, Sobs_list):
    """
    The log likelihood given the data

    """

    # Construct covariance
    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    sigma = M * M.transpose()

    # Compute the loglikelihood
    lnL = 0
    for i in range(len(s2_list)):

        # Set stuff
        s2 = s2_list[i]
        xbar = xbar_list[i]
        ctot = ctot_list[i]
        Sobs = Sobs_list[i]

        # Get the change of basis operation
        R = compute_change_of_basis_operation(s0, s2)

        # Compute rotated sigma
        S = R * sigma * R.transpose()
        S11 = matrix.sqr((S[0], S[1], S[3], S[4]))
        S12 = matrix.col((S[2], S[5]))
        S21 = matrix.col((S[6], S[7])).transpose()
        S22 = S[8]

        # Compute rotated mu
        mu = R * s2
        assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]

        # Compute conditional mean and covariance
        mubar = mu1 + S12 * (1 / S22) * (s0.length() - mu2)
        Sbar = S11 - S12 * (1 / S22) * S21

        # Compute the log likelihood
        try:
            Sobs = matrix.sqr(Sobs)
            B1 = log(S22)
            B2 = (1 / S22) * (s0.length() - mu2) ** 2
            A1 = log(Sbar.determinant()) * ctot
            A2 = ((Sbar.inverse()) * (ctot * Sobs)).trace()
            A3 = (
                (Sbar.inverse()) * (ctot * (xbar - mubar) * (xbar - mubar).transpose())
            ).trace()
            lnL += -0.5 * (A1 + A2 + A3 + (B1 + B2))
        except Exception:
            raise
            lnL += -1e15

    # print tuple(sigma), lnL
    return lnL


class Target(object):
    def __init__(self, s0, s2_list, xbar_list, ctot_list, Sobs_list):
        self.s0 = s0
        self.s2_list = s2_list
        self.xbar_list = xbar_list
        self.ctot_list = ctot_list
        self.Sobs_list = Sobs_list

    def target(self, params):
        lnL = log_likelihood(
            params,
            self.s0,
            self.s2_list,
            self.xbar_list,
            self.ctot_list,
            self.Sobs_list,
        )
        score = -lnL
        return score


def test_ideal():

    numpy.random.seed(100)

    # The beam vector
    s0 = matrix.col((0, 0, 1))

    # The covariance matrix
    sigma = matrix.sqr((1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 3e-6))

    # The number of reflections
    N = 100

    # Generate a load of reflections
    s2_list, ctot_list, xbar_list, Sobs_list = generate_simple(s0, sigma, N=N)

    # Starting values for simplex
    values = flex.double((sqrt(1e-6), 0, sqrt(1e-6), 0, 0, sqrt(1e-6)))
    offset = flex.double([sqrt(1e-7) for v in values])

    # Do the simplex optimization
    optimizer = SimpleSimplex(
        values, offset, Target(s0, s2_list, xbar_list, ctot_list, Sobs_list), 2000
    )
    params = optimizer.get_solution()

    # Create the covariance matrix
    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    sigma = M * M.transpose()

    print(sigma)

    expected = matrix.sqr(
        (
            9.91047018199e-07,
            -1.98078253593e-09,
            2.27093231797e-09,
            -1.98078253593e-09,
            1.98335548957e-06,
            1.88051940862e-08,
            2.27093231797e-09,
            1.88051940862e-08,
            2.99885951955e-06,
        )
    )
    assert all(1e6 * abs(a - b) < 1e-7 for a, b in zip(sigma, expected))

    print("OK")


def test_binned():

    numpy.random.seed(100)

    # The beam vector
    s0 = matrix.col((0, 0, 1))

    # The covariance matrix
    sigma = matrix.sqr((1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 3e-6))

    # The number of reflections
    N = 100

    # Generate a load of reflections
    s2_list, ctot_list, xbar_list, Sobs_list = generate_simple_binned(s0, sigma, N=N)

    # Starting values for simplex
    values = flex.double((sqrt(1e-6), 0, sqrt(1e-6), 0, 0, sqrt(1e-6)))
    offset = flex.double([sqrt(1e-7) for v in values])

    # Do the simplex optimization
    optimizer = SimpleSimplex(
        values, offset, Target(s0, s2_list, xbar_list, ctot_list, Sobs_list), 2000
    )
    params = optimizer.get_solution()

    # Create the covariance matrix
    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    sigma = M * M.transpose()

    # print sigma

    expected = matrix.sqr(
        (
            1.04779078869e-06,
            -1.38584634956e-09,
            -5.72912737193e-09,
            -1.38584634956e-09,
            2.09196281286e-06,
            3.30231414058e-08,
            -5.72912737193e-09,
            3.30231414058e-08,
            3.1533805736e-06,
        )
    )

    assert all(1e6 * abs(a - b) < 1e-7 for a, b in zip(sigma, expected))

    print("OK")


@pytest.mark.xfail(
    reason="generate_simple does not generate everything needed - outdated code?"
)
def test_ml_target_class():
    class SimplexTarget(object):
        def __init__(self, s0, s2_list, ctot_list, xbar_list, Sobs_list):
            self.s0 = s0
            self.s2_list = s2_list
            self.xbar_list = xbar_list
            self.ctot_list = ctot_list
            self.Sobs_list = Sobs_list

        def target(self, params):

            parameterisation = Simple6MosaicityParameterisation(params)
            # model, s0, sp_list, h_list, ctot_list, mobs_list, sobs_list
            t = MaximumLikelihoodTarget(
                parameterisation,
                self.s0,
                self.s2_list,
                self.ctot_list,
                self.xbar_list,
                self.Sobs_list,
            )

            lnL = t.log_likelihood()

            # print tuple(parameterisation.sigma()), lnL

            return -lnL

    numpy.random.seed(100)

    # The beam vector
    s0 = matrix.col((0, 0, 1))

    # The covariance matrix
    sigma = matrix.sqr((1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 3e-6))

    # The number of reflections
    N = 100

    # Generate a load of reflections
    s2_list, ctot_list, xbar_list, Sobs_list = generate_simple(s0, sigma, N=N)

    Sobs_list = flex.double(Sobs_list)

    # Starting values for simplex
    values = flex.double((sqrt(1e-6), 0, sqrt(1e-6), 0, 0, sqrt(1e-6)))
    offset = flex.double([sqrt(1e-7) for v in values])

    # Do the simplex optimization
    optimizer = SimpleSimplex(
        values,
        offset,
        SimplexTarget(s0, s2_list, ctot_list, xbar_list, Sobs_list),
        2000,
    )
    params = optimizer.get_solution()

    # Create the covariance matrix
    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    sigma = M * M.transpose()

    expected = matrix.sqr(
        (
            9.91047018199e-07,
            -1.98078253593e-09,
            2.27093231797e-09,
            -1.98078253593e-09,
            1.98335548957e-06,
            1.88051940862e-08,
            2.27093231797e-09,
            1.88051940862e-08,
            2.99885951955e-06,
        )
    )
    assert all(1e6 * abs(a - b) < 1e-7 for a, b in zip(sigma, expected))

    print("OK")


@pytest.mark.xfail(
    reason="generate_simple does not generate everything needed - outdated code?"
)
def test_ml_target_class_2():

    numpy.random.seed(100)

    # The beam vector
    s0 = matrix.col((0, 0, 1))

    # The covariance matrix
    sigma = matrix.sqr((1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 3e-6))

    # The number of reflections
    N = 100

    # Generate a load of reflections
    s2_list, ctot_list, xbar_list, Sobs_list = generate_simple(s0, sigma, N=N)

    Sobs_list = flex.double(Sobs_list)

    parameterisation = Simple6MosaicityParameterisation((1, 0, 1, 0, 0, 1))

    data = ProfileRefinerData(s0, s2_list, ctot_list, xbar_list, Sobs_list)
    refiner = ProfileRefiner(parameterisation, data)
    refiner.refine()
    params = refiner.parameters

    # Create the covariance matrix
    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    sigma = M * M.transpose()

    # print sigma

    expected = matrix.sqr(
        (
            9.91048657253e-07,
            -1.9828296735e-09,
            2.25787032072e-09,
            -1.9828296735e-09,
            1.98334108426e-06,
            1.88097904832e-08,
            2.25787032072e-09,
            1.88097904832e-08,
            2.99884748097e-06,
        )
    )

    assert all(1e6 * abs(a - b) < 1e-7 for a, b in zip(sigma, expected))

    print("OK")


if __name__ == "__main__":
    test_ideal()
    test_binned()
    test_ml_target_class()
    test_ml_target_class_2()
