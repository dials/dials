from __future__ import division, print_function

from math import log
from os.path import join

import numpy.random
from numpy.random import choice as sample

from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix

from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.potato.parameterisation import (
    ModelState,
    Simple6MosaicityParameterisation,
)
from dials.algorithms.profile_model.potato.refiner import Refiner as ProfileRefiner
from dials.algorithms.profile_model.potato.refiner import (
    RefinerData as ProfileRefinerData,
)
from dials.algorithms.profile_model.potato.util.generate_simple import (
    generate_from_reflections,
    generate_from_reflections_binned,
)
from dials.array_family import flex


def log_likelihood(params, s0, s2_list, xbar_list, ctot_list, Sobs_list, test=0):

    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )

    sigma = M * M.transpose()

    lnL = 0
    for i in range(len(s2_list)):

        s2 = s2_list[i]
        xbar = xbar_list[i]
        ctot = ctot_list[i]
        Sobs = Sobs_list[i]

        R = compute_change_of_basis_operation(s0, s2)
        S = R * sigma * R.transpose()
        S11 = matrix.sqr((S[0], S[1], S[3], S[4]))
        S12 = matrix.col((S[2], S[5]))
        S21 = matrix.col((S[6], S[7])).transpose()
        S22 = S[8]
        mu = R * s2
        assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]

        mubar = mu1 + S12 * (1 / S22) * (s0.length() - mu2)
        Sbar = S11 - S12 * (1 / S22) * S21
        try:
            Sobs = matrix.sqr(Sobs)
            B1 = log(S22)
            B2 = (1 / S22) * (s0.length() - mu2) ** 2
            A1 = log(Sbar.determinant()) * ctot
            A2 = ((Sbar.inverse()) * (Sobs)).trace()
            A3 = (
                (Sbar.inverse()) * (ctot * (xbar - mubar) * (xbar - mubar).transpose())
            ).trace()
            if test:
                lnL += -0.5 * (A1 + A2 + A3 + ctot * (B1 + B2))
            else:
                lnL += -0.5 * (A1 + A2 + A3 + (B1 + B2))
        except Exception:
            lnL += -1e15

    print(tuple(sigma), lnL)
    return lnL


class Target(object):
    def __init__(self, s0, s2_list, xbar_list, ctot_list, Sobs_list, test=1):
        self.s0 = s0
        self.s2_list = s2_list
        self.xbar_list = xbar_list
        self.ctot_list = ctot_list
        self.Sobs_list = Sobs_list
        self.test = test

    def target(self, params):

        lnL = log_likelihood(
            params,
            self.s0,
            self.s2_list,
            self.xbar_list,
            self.ctot_list,
            self.Sobs_list,
            test=self.test,
        )

        score = -lnL
        return score


def generate_observations2(experiments, reflections, sigma):

    A = matrix.sqr(experiments[0].crystal.get_A())
    s0 = matrix.col(experiments[0].beam.get_s0())

    s2_obs = flex.vec3_double()
    for i in range(len(reflections)):

        h = matrix.col(reflections[i]["miller_index"])

        r = A * h
        s2 = s0 + r

        s2_obs.append(s2)

    reflections["s2"] = s2_obs
    return reflections


def test_ideal(dials_regression):

    numpy.random.seed(100)

    # Ensure we have a data block
    # experiments = ExperimentListFactory.from_json_file("experiments.json")
    # was being loaded locally, assumed to be file moved into dials_regression
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 1.0), deg=True)
    experiments[0].beam.set_s0((0, 0, -1))

    s0 = matrix.col(experiments[0].beam.get_s0())

    # The predicted reflections
    reflections = flex.reflection_table.from_predictions_multi(experiments, padding=4)
    print(len(reflections))

    sigma = matrix.sqr((1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 3e-6))

    reflections = generate_observations2(experiments, reflections, sigma)

    s2_list, h_list, ctot_list, xbar_list, Sobs_list = generate_from_reflections(
        s0, sigma, reflections
    )

    index = sample(range(len(s2_list)), 200)

    def select_sample(d, index):
        return [d[i] for i in index]

    s2_list = select_sample(s2_list, index)
    h_list = select_sample(h_list, index)
    ctot_list = select_sample(ctot_list, index)
    xbar_list = select_sample(xbar_list, index)
    Sobs_list = select_sample(Sobs_list, index)

    print("Using %d reflections: " % len(s2_list))

    # values = flex.double((sqrt(1.1e-6), 0, sqrt(2.1e-6), 0, 0, sqrt(3.1e-6)))
    # offset = flex.double([sqrt(1e-7) for v in values])

    parameterisation = Simple6MosaicityParameterisation(flex.double([1, 0, 1, 0, 0, 1]))
    Sobs_list = flex.double(Sobs_list)
    data = ProfileRefinerData(s0, s2_list, h_list, ctot_list, xbar_list, Sobs_list)
    state = ModelState(experiments[0], parameterisation)
    refiner = ProfileRefiner(state, data)
    refiner.refine()
    params = refiner.state.M_parameterisation.parameters
    # optimizer = SimpleSimplex(
    #   values,
    #   offset,
    #   Target(
    #     s0,
    #     s2_list,
    #     xbar_list,
    #     ctot_list,
    #     Sobs_list,
    #     test=0), 2000)
    # params = optimizer.get_solution()

    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )

    sigma = M * M.transpose()
    print(sigma)

    expected = matrix.sqr(
        (
            1.0030467686e-06,
            -1.98473936999e-09,
            -8.60673302905e-10,
            -1.98473936999e-09,
            2.00630994244e-06,
            -1.64963854836e-08,
            -8.60673302905e-10,
            -1.64963854836e-08,
            2.97450815302e-06,
        )
    )

    assert all(1e6 * abs(a - b) < 1e-7 for a, b in zip(sigma, expected))

    print("OK")


def test_binned(dials_regression):

    numpy.random.seed(100)

    # Ensure we have a data block
    # experiments = ExperimentListFactory.from_json_file("experiments.json")
    # was being loaded locally, assumed to be file moved into dials_regression
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 1.0), deg=True)
    experiments[0].beam.set_s0((0, 0, -1))

    s0 = matrix.col(experiments[0].beam.get_s0())

    # The predicted reflections
    reflections = flex.reflection_table.from_predictions_multi(experiments, padding=4)
    print(len(reflections))

    sigma = matrix.sqr((1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 3e-6))

    reflections = generate_observations2(experiments, reflections, sigma)

    s2_list, hlist, ctot_list, xbar_list, Sobs_list = generate_from_reflections_binned(
        s0, sigma, reflections
    )

    index = sample(range(len(s2_list)), 200)

    def select_sample(d, index):
        return [d[i] for i in index]

    s2_list = select_sample(s2_list, index)
    hlist = select_sample(hlist, index)
    ctot_list = select_sample(ctot_list, index)
    xbar_list = select_sample(xbar_list, index)
    Sobs_list = select_sample(Sobs_list, index)

    print("Using %d reflections: " % len(s2_list))

    # values = flex.double((sqrt(1e-6), 0, sqrt(2e-6), 0, 0, sqrt(3e-6)))
    # offset = flex.double([sqrt(1e-7) for v in values])

    parameterisation = Simple6MosaicityParameterisation(flex.double([1, 0, 1, 0, 0, 1]))
    Sobs_list = flex.double(Sobs_list)
    data = ProfileRefinerData(s0, s2_list, hlist, ctot_list, xbar_list, Sobs_list)
    state = ModelState(
        experiments[0], parameterisation, fix_unit_cell=True, fix_orientation=True
    )
    refiner = ProfileRefiner(state, data)
    refiner.refine()
    params = refiner.state.M_parameterisation.parameters
    # optimizer = SimpleSimplex(
    #   values,
    #   offset,
    #   Target(
    #     s0,
    #     s2_list,
    #     xbar_list,
    #     ctot_list,
    #     Sobs_list,
    #     test=0), 2000)
    # params = optimizer.get_solution()

    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )

    sigma = M * M.transpose()
    print(sigma)

    expected = matrix.sqr(
        (
            1.07025484551e-06,
            1.30518861783e-09,
            -1.72635922351e-09,
            1.30518861783e-09,
            2.10252906788e-06,
            -1.64646310672e-08,
            -1.72635922351e-09,
            -1.64646310672e-08,
            3.12149393966e-06,
        )
    )
    assert all(1e6 * abs(a - b) < 1e-7 for a, b in zip(sigma, expected))

    print("OK")


if __name__ == "__main__":
    test_ideal()
    test_binned()
