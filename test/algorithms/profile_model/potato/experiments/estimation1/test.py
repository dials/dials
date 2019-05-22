from __future__ import division, print_function

from collections import defaultdict
from math import cos, log, pi, sin
from random import sample, uniform

from scitbx import matrix

from dials.array_family import flex
from dials_scratch.jmp.potato.parameterisation import SimpleMosaicityParameterisation
from dials_scratch.jmp.potato.profile_refiner import ProfileRefiner, ProfileRefinerData
from dials_scratch.jmp.potato.util.generate_simple import generate_simple


def generate_sigma():
    """
    Generate a random sigma

    """

    # Random eigenvalues
    L = matrix.diag((uniform(1e-5, 1e-4), uniform(1e-5, 1e-4), uniform(1e-5, 1e-4)))

    # Random angles
    ax = uniform(0, 2 * pi)
    ay = uniform(0, 2 * pi)
    az = uniform(0, 2 * pi)

    # X rotation
    Rx = matrix.sqr((1, 0, 0, 0, cos(ax), -sin(ax), 0, sin(ax), cos(ax)))

    # Y rotation
    Ry = matrix.sqr((cos(ay), 0, sin(ay), 0, 1, 0, -sin(ay), 0, cos(ay)))

    # Z rotation
    Rz = matrix.sqr((cos(az), -sin(az), 0, sin(az), cos(az), 0, 0, 0, 1))

    # Eigen vectors
    Q = Rz * Ry * Rx

    # Compute sigma
    sigma = Q * L * Q.inverse()

    return sigma


def kl_divergence(A, B):
    return 0.5 * (
        (B.inverse() * A).trace() - 3 + log(B.determinant() / A.determinant())
    )


def run():
    # np.random.seed(1)
    # seed(1)

    # The number of times to run each experiment
    m = 10

    # The results
    results = defaultdict(list)

    # Generate a covariance matrix
    print("Generating covariance matrix")
    sigma = generate_sigma()

    # Repeat the experiment m times
    for i in range(m):

        # Set the beam vector
        s0 = matrix.col((0, 0, 1))
        # sigma = matrix.sqr((2.7e-05, -2.33e-05, -1.26e-06,
        #                     -2.33e-05, 6.08e-05, -1.01e-05,
        #                     -1.26e-06, -1.01e-05, 4.41e-05))
        print("Known Sigma")
        print(
            "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, )" % tuple(sigma)
        )

        # Generate some samples
        print("Generating samples")
        samples = list(zip(*generate_simple(s0, sigma, 1000)))

        # Iterate over the number of samples to use (logarithmically between 10 and
        # 1000).
        # A = log(10)
        # B = log(1001)
        # N = (B - A) / 10.0
        # for t in np.arange(A, log(1001), N):
        #   n = int(floor(exp(t)))
        # N_list = [900]
        N_list = list(range(3, 20))
        N_list = (
            list(range(3, 20)) + list(range(20, 200, 10)) + list(range(200, 1000, 100))
        )
        for n in N_list:

            # Select a sample
            print("Selecting %d samples" % n)
            subsample = sample(samples, n)

            # Estimate the parameters
            print("Estimating parameters")
            s2_list, ctot_list, mobs_list, Sobs_list = zip(*subsample)
            Sobs_list = flex.double(Sobs_list)
            refiner = ProfileRefiner(
                SimpleMosaicityParameterisation((1, 0, 1, 0, 0, 1)),
                ProfileRefinerData(s0, s2_list, ctot_list, mobs_list, Sobs_list),
            )
            refiner.refine()
            params = refiner.parameters

            # Compute sigma from the parameters
            M = matrix.sqr(
                (
                    params[0],
                    0,
                    0,
                    params[1],
                    params[2],
                    0,
                    params[3],
                    params[4],
                    params[5],
                )
            )
            sigma_cal = M * M.transpose()

            # Compute the KL divergence
            kl = kl_divergence(sigma, sigma_cal)

            print("Calculated Sigma")
            print(
                "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, )"
                % tuple(sigma_cal)
            )
            print(n, kl)
            print("")

            results[n].append(kl)

        print("Known Sigma")
        print(
            "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, )" % tuple(sigma)
        )

    with open("result.txt", "w") as outfile:
        for n in sorted(results.keys()):
            print("%d %.6f" % (n, sum(results[n]) / len(results[n])), file=outfile)


if __name__ == "__main__":
    run()
