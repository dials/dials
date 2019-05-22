import pytest
from scitbx import matrix
from dials.array_family import flex

# from matplotlib import pylab
from math import sqrt, exp
from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)


def test_ewald_spread_approximation():

    # The wavelength spread
    sigma_l = 0.001
    sigma_r = matrix.sqr((1e-7, 0, 0, 0, 1e-7, 0, 0, 0, 5e-8))

    # The beam vector
    s0 = matrix.col((0, 0, -1))

    # A test point
    s1 = matrix.col((0, 1, -1)).normalize() * s0.length()
    r0 = s1 - s0

    R = compute_change_of_basis_operation(s0, s1)

    # A 3D box around the test point
    d = 0.001
    box = (s1[0] - d, s1[0] + d, s1[1] - d, s1[1] + d, s1[2] - d, s1[2] + d)

    grid_size = 51
    step_size = (2 * d) / float(grid_size)

    expected = flex.double(flex.grid(grid_size, grid_size, grid_size))
    approximated = flex.double(flex.grid(grid_size, grid_size, grid_size))
    approximated_dist = flex.double(flex.grid(grid_size, grid_size, grid_size))
    rlp_dist = flex.double(flex.grid(grid_size, grid_size, grid_size))
    approximated_product_dist = flex.double(flex.grid(grid_size, grid_size, grid_size))

    for k in range(grid_size):
        for j in range(grid_size):
            for i in range(grid_size):
                s = matrix.col(
                    (
                        box[0] + step_size * i,
                        box[2] + step_size * j,
                        box[4] + step_size * k,
                    )
                )
                r = s - s0
                expected[k, j, i] = (-2 / s0.length()) * (s0.dot(r)) / (r.dot(r))

                s_e = s.normalize() * s0.length()
                r_e = s_e - s0
                z = ((r + s0) * s0.length()).dot((s0 + r_e))

                approximated[k, j, i] = (1.0 / s0.length()) - 2 * (z - s0.length()) / (
                    r_e.dot(r_e)
                )

                sigma_E2 = sigma_l ** 2 * (r_e.dot(r_e) / 2.0) ** 2
                sigma_E = sqrt(sigma_E2)

                approximated_dist[k, j, i] = exp(
                    -0.5 * (z - s0.length()) ** 2 / sigma_E ** 2
                )

                rlp_dist[k, j, i] = exp(
                    -0.5 * ((r - r0).transpose() * sigma_r.inverse() * (r - r0))[0]
                )

    assert max(flex.abs((expected - approximated) / expected)) < 1e-5

    # pylab.imshow(expected.as_numpy_array()[:,:,25])
    # # # pylab.imshow(expected.as_numpy_array()[:,25,:])
    # # # # pylab.imshow(expected.as_numpy_array()[25,:,:])
    # pylab.show()
    # pylab.imshow(approximated.as_numpy_array()[:,:,25])
    # # # pylab.imshow(approximated.as_numpy_array()[:,25,:])
    # # # # pylab.imshow(approximated.as_numpy_array()[25,:,:])
    # pylab.show()
    # pylab.imshow((expected-approximated).as_numpy_array()[:,:,25])
    # # # pylab.imshow(approximated.as_numpy_array()[:,25,:])
    # # # # pylab.imshow(approximated.as_numpy_array()[25,:,:])
    # pylab.show()

    expected_dist = flex.exp(-0.5 * (expected - 1.0) ** 2 / sigma_l ** 2)
    # pylab.imshow(expected_dist.as_numpy_array()[:,:,25])
    # pylab.show()
    # pylab.imshow(approximated_dist.as_numpy_array()[:,:,25])
    # pylab.show()
    # pylab.imshow(((expected_dist - approximated_dist)/expected_dist).as_numpy_array()[:,:,25])
    # pylab.show()
    assert max(flex.abs((expected_dist - approximated_dist) / expected_dist)) < 0.02

    expected_product_dist = expected_dist * rlp_dist

    # Set SigmaE once
    r_e = s1 - s0
    sigma_E2 = sigma_l ** 2 * (r_e.dot(r_e) / 2.0) ** 2
    mu = R * s1
    Sigma = R * sigma_r * R.transpose()
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]
    Sigma11 = matrix.sqr((Sigma[0], Sigma[1], Sigma[3], Sigma[4]))
    Sigma12 = matrix.col((Sigma[2], Sigma[5]))
    Sigma21 = matrix.col((Sigma[6], Sigma[7])).transpose()
    Sigma22 = Sigma[8]
    kappa = sigma_E2 / (Sigma22 + sigma_E2)
    p2 = (mu2 * sigma_E2 + s0.length() * Sigma22) / (Sigma22 + sigma_E2)
    p1 = mu1 + Sigma12 / Sigma22 * (p2 - mu2)
    p = matrix.col((p1[0], p1[1], p2))
    P11 = Sigma11 - Sigma12 * Sigma21 * (1 - kappa) / Sigma22
    P12 = kappa * Sigma12
    P21 = kappa * Sigma21
    P22 = kappa * Sigma22
    P = matrix.sqr(
        (P11[0], P11[1], P12[0], P11[2], P11[3], P12[1], P21[0], P21[1], P22)
    )

    # Compute the approximated produce dist
    for k in range(grid_size):
        for j in range(grid_size):
            for i in range(grid_size):
                s = matrix.col(
                    (
                        box[0] + step_size * i,
                        box[2] + step_size * j,
                        box[4] + step_size * k,
                    )
                )
                r = s - s0
                x = R * s
                approximated_product_dist[k, j, i] = exp(
                    -0.5 * ((x - p).transpose() * P.inverse() * (x - p))[0]
                )

    # pylab.imshow(expected_product_dist.as_numpy_array()[:,:,25])
    # pylab.show()
    # pylab.imshow(approximated_product_dist.as_numpy_array()[:,:,25])
    # pylab.show()
    # pylab.imshow(((expected_product_dist - approximated_product_dist)/expected_product_dist).as_numpy_array()[:,:,25])
    # pylab.show()

    assert (
        max(
            flex.abs(
                (approximated_product_dist - expected_product_dist)
                / expected_product_dist
            )
        )
    ) < 0.02

    # Compute the covariance matrix
    expected_covariance = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
    approximated_covariance = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
    n1 = flex.sum(expected_product_dist)
    n2 = flex.sum(approximated_product_dist)
    for k in range(grid_size):
        for j in range(grid_size):
            for i in range(grid_size):
                s = matrix.col(
                    (
                        box[0] + step_size * i,
                        box[2] + step_size * j,
                        box[4] + step_size * k,
                    )
                )
                r = s - s0
                x = R * s
                expected_covariance += (
                    expected_product_dist[k, j, i] * (x - p) * (x - p).transpose()
                )
                approximated_covariance += (
                    approximated_product_dist[k, j, i] * (x - p) * (x - p).transpose()
                )
    expected_covariance /= n1
    approximated_covariance /= n2

    assert tuple(expected_covariance) == pytest.approx(tuple(approximated_covariance))
    assert tuple(approximated_covariance) == pytest.approx(tuple(P), rel=0.1)

    expected_product_dist /= flex.sum(expected_product_dist)
    approximated_product_dist /= flex.sum(approximated_product_dist)

    # Compute the distribution of diffracted beam vectors
    identity = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
    rc = p
    J = (
        identity
        - 0.5
        * s0
        * ((2 * (rc.dot(s0)) * rc - (rc.dot(rc)) * s0) / (rc.dot(s0)) ** 2).transpose()
    )
    sc = -0.5 * (rc.dot(rc) / s0.dot(rc)) * s0 + rc
    q = sc
    Q = J * P * J.transpose()
    Q11 = matrix.sqr((Q[0], Q[1], Q[3], Q[4]))
    # q1 = matrix.col((q[0], q[1]))

    # diffracted_dist = flex.double(flex.grid(grid_size//2, grid_size//2))
    # diffracted_dist2 = flex.double(flex.grid(grid_size//4, grid_size//4, grid_size//4))

    covariance_diffracted_dist = matrix.sqr((0, 0, 0, 0))
    for k in range(grid_size):
        for j in range(grid_size):
            for i in range(grid_size):
                w1 = expected_product_dist[k, j, i]
                s = matrix.col(
                    (
                        box[0] + step_size * i,
                        box[2] + step_size * j,
                        box[4] + step_size * k,
                    )
                )
                r = s - s0
                l = (-2 / s0.length()) * (s0.dot(r)) / (r.dot(r))
                s00 = s0.normalize() / l
                s11 = s00 + r

                # kk = int((s11[2] - (box[4]-d)) / (8*step_size))
                # jj = int((s11[1] - (box[2]-d)) / (8*step_size))
                # ii = int((s11[0] - (box[0]-d)) / (8*step_size))

                # if all(item >= 0 and item < grid_size//4 for item in [kk, jj, ii]):
                #     diffracted_dist2[kk,jj,ii] += w1

                assert s11.length() == pytest.approx(s00.length())
                x = R * s  # 11
                z = matrix.col((0, 0, 1.0))
                t = 1.0 / x.dot(z)
                xx = t * x
                q = matrix.col((xx[0], xx[1]))
                covariance_diffracted_dist += w1 * q * q.transpose()
                # ii = int((xx[0]+0.001)/(2*step_size))
                # jj = int((xx[1]+0.001)/(2*step_size))
                # if ii >= 0 and jj >= 0 and ii < grid_size//2 and jj < grid_size//2:
                #     diffracted_dist[jj,ii] += w1
    covariance_diffracted_dist /= flex.sum(expected_product_dist)

    assert tuple(covariance_diffracted_dist) == pytest.approx(tuple(Q11), rel=0.1)

    # print("")
    # print(covariance_diffracted_dist)
    # print(Q11)

    # print(var)
    # from matplotlib import pylab
    # pylab.imshow(diffracted_dist2.as_numpy_array()[:,:,12])
    # pylab.show()
