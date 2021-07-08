from __future__ import division, print_function

from math import sqrt
from os.path import join
from random import sample, seed

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix

from dials.algorithms.profile_model.potato.model import (
    Simple6ProfileModel,
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.potato.refiner import Refiner
from dials.array_family import flex


def compute_mean_plane(mu, sigma, s0):
    z = matrix.col((0, 0, 1))
    R = compute_change_of_basis_operation(s0, mu)

    sigma_1 = R * sigma * R.transpose()
    mu_1 = R * mu

    assert abs(1 - mu_1.normalize().dot(z)) < 1e-7

    # sigma_11 = matrix.sqr((sigma_1[0], sigma_1[1], sigma_1[3], sigma_1[4]))
    sigma_12 = matrix.col((sigma_1[2], sigma_1[5]))
    # sigma_21 = matrix.col((sigma_1[3], sigma_1[7])).transpose()
    sigma_22 = sigma_1[8]
    mu1 = matrix.col((mu_1[0], mu_1[1]))
    mu2 = mu_1[2]
    mu_new_1 = mu1 + sigma_12 * (1 / sigma_22) * (s0.length() - mu2)
    v = matrix.col((mu_new_1[0], mu_new_1[1], s0.length())).normalize() * s0.length()
    x_new = R.transpose() * v
    return x_new


def generate_observations(experiments, reflections, sigma):

    A = matrix.sqr(experiments[0].crystal.get_A())
    s0 = matrix.col(experiments[0].beam.get_s0())

    s1_obs = flex.vec3_double()
    s2_obs = flex.vec3_double()
    for i in range(len(reflections)):

        h = matrix.col(reflections[i]["miller_index"])

        r = A * h
        s2 = s0 + r

        s1 = compute_mean_plane(s2, sigma, s0)

        s1_obs.append(s1)
        s2_obs.append(s2)

    return s1_obs, s2_obs


@pytest.mark.xfail(reason="outdated code")
def test_simplex(dials_regression):

    from dials.array_family import flex

    seed(0)

    # Ensure we have a data block
    # experiments = ExperimentListFactory.from_json_file("experiments.json")
    # was loading a local file, assumed to be that moved to dials_regression
    experiments = ExperimentListFactory.from_json_file(
        join(dials_regression, "potato_test_data", "experiments.json")
    )
    experiments[0].scan.set_oscillation((0, 1), deg=True)
    # experiments[0].scan = experiments[0].scan[0:1]
    # experiments[0].imageset = experiments[0].imageset[0:1]

    # The predicted reflections
    reflections = flex.reflection_table.from_predictions_multi(experiments, padding=1)

    # Select only those within 1 deg
    x, y, z = reflections["xyzcal.px"].parts()
    selection = flex.abs(z) < 1
    reflections = reflections.select(selection)

    selection = flex.size_t(sample(range(len(reflections)), 1000))
    reflections = reflections.select(selection)

    parameters = (sqrt(0.0001), 0, sqrt(0.0002), 0, 0, sqrt(0.0003))
    M = matrix.sqr(
        (
            parameters[0],
            0,
            0,
            parameters[1],
            parameters[2],
            0,
            parameters[3],
            parameters[4],
            parameters[5],
        )
    )
    sigma = M * M.transpose()

    print(sigma)

    # Generate observed positions
    s1_obs, s2_obs = generate_observations(experiments, reflections, sigma)

    angles = []
    for s1, s2 in zip(s1_obs, s2_obs):
        a = matrix.col(s1).angle(matrix.col(s2), deg=True)
        angles.append(a)
    print("Mean angle between s1 and s2 %f degrees " % (sum(angles) / len(angles)))

    # Do the ray intersection
    reflections["s1_obs"] = s1_obs
    reflections["s1"] = s2_obs

    xyzobs = flex.vec3_double()
    xyzobspx = flex.vec3_double()
    for j in range(len(s1_obs)):
        mm = experiments[0].detector[0].get_ray_intersection(s1_obs[j])
        px = experiments[0].detector[0].millimeter_to_pixel(mm)
        xyzobs.append((mm[0], mm[1], 0))
        xyzobspx.append((px[0], px[1], 0))

    reflections["xyzobs.mm.value"] = xyzobs
    reflections["xyzobs.px.value"] = xyzobspx

    # Offset the crystal orientation matrix
    U = matrix.sqr(experiments[0].crystal.get_U())

    print(
        "Original orientation: ",
        "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f)" % tuple(U),
    )

    m2 = matrix.col(experiments[0].goniometer.get_rotation_axis())
    R = m2.axis_and_angle_as_r3_rotation_matrix(angle=0.5, deg=True)
    experiments[0].crystal.set_U(R * U)
    model = Simple6ProfileModel.from_sigma(sigma)
    # model = Simple6MosaicityModel(sigma) - old, assume above is what updated
    # version is

    # Do the refinement
    from unittest.mock import Mock

    params = Mock()
    params.refinement.n_cycles = 3
    params.profile.rlp_mosaicity.model = "simple6"
    refiner = Refiner(experiments, reflections, model, params)

    crystal = refiner.experiments.crystal

    U_old = U
    U = crystal.get_U()

    print(
        "Refined orientation: ",
        "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f)" % tuple(U),
    )

    assert all(abs(u1 - u2) < 1e-7 for u1, u2 in zip(U_old, U))

    print("OK")


if __name__ == "__main__":
    test_simplex()
