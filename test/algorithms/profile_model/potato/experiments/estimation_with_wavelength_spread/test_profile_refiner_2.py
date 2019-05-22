from __future__ import division, print_function

from math import log, sqrt

import numpy.random
from scipy.optimize import minimize

from scitbx import matrix

from dials.array_family import flex
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from dials_scratch.jmp.potato.util.generate_simple import (
    generate_with_wavelength_spread2,
)


def compute_beam_vector_rotation(s0):
    """
    Construct a matrix to rotate whole coordinate system so that beam is along z

    """
    z_axis = s0.normalize()
    n = [i for i in range(3) if z_axis[i] != 0][-1]
    temp = [0, 0, 0]
    temp[(n + 1) % 3] = 1
    y_axis = z_axis.cross(matrix.col(temp)).normalize()
    x_axis = y_axis.cross(z_axis).normalize()
    R = matrix.sqr(x_axis.elems + y_axis.elems + z_axis.elems)
    assert abs((R * s0.normalize()).dot(matrix.col((0, 0, 1))) - 1) < 1e-7
    return R


def matrix_from_params(params):
    # L = matrix.diag((
    #   (0.001*atan(params[0]) / (pi / 2))**2,
    #   (0.001*atan(params[1]) / (pi / 2))**2,
    #   (0.001*atan(params[2]) / (pi / 2))**2))
    # Rx = matrix.sqr((
    #   1, 0, 0,
    #   0, cos(params[3]), -sin(params[3]),
    #   0, sin(params[3]), cos(params[3])))
    # Ry = matrix.sqr((
    #   cos(params[4]), 0, sin(params[4]),
    #   0, 1, 0,
    #   -sin(params[4]), 0, cos(params[4])))
    # Rz = matrix.sqr((
    #   cos(params[5]), -sin(params[5]), 0,
    #   sin(params[5]), cos(params[5]), 0,
    #   0, 0, 1))
    # Q = (Rz*Ry*Rx).transpose()
    # return Q*L*Q.transpose()
    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    return M * M.transpose()


def log_likelihood(params, s0, s2_list, xbar_list, ctot_list, Sobs_list):
    """
    The log likelihood given the data

    """

    # Construct covariance
    spot_covariance = matrix_from_params(params)
    # M = matrix.sqr((
    #   params[0], 0, 0,
    #   params[1], params[2], 0,
    #   params[3], params[4], params[5]))
    # spot_covariance = M*M.transpose()
    # spot_covariance = matrix.sqr((
    #   1e-7, 0, 0,
    #   0, 2e-7, 0,
    #   0, 0, 3e-7))
    # Construct the wavelength spread
    # wavelength_variance = params[6]**2

    wavelength_variance = params[6] ** 2  # (0.05*atan(params[6]) / (pi / 2))**2
    # wavelength_variance = (0.05*atan(params[0]) / (pi / 2))**2

    R0 = compute_beam_vector_rotation(s0)

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
        S = R * spot_covariance * R.transpose()

        # Partition matrix
        S11 = matrix.sqr((S[0], S[1], S[3], S[4]))
        S12 = matrix.col((S[2], S[5]))
        S21 = matrix.col((S[6], S[7])).transpose()
        S22 = S[8]

        # Compute the local spread in energy
        q0 = s2.normalize() * s0.length() - s0
        wavelength_variance_local = wavelength_variance * (q0.dot(q0) / 2) ** 2

        # Apply the wavelength spread to the sigma
        Sp_22 = S22 * wavelength_variance_local / (S22 + wavelength_variance_local)
        Sp_12 = S12 * Sp_22 / S22
        Sp_11 = S11 - (S12 * (1 / S22) * S21) * (1 - Sp_22 / S22)

        Sigma3 = matrix.sqr(
            (
                Sp_11[0],
                Sp_11[1],
                Sp_12[0],
                Sp_11[2],
                Sp_11[3],
                Sp_12[1],
                Sp_12[0],
                Sp_12[1],
                Sp_22,
            )
        )

        mu = R * s2
        assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7

        # Apply the wavelength spread to the mean
        mu1 = mu
        mu2 = matrix.col((0, 0, s0.length()))
        z0 = (mu1[2] * wavelength_variance_local + mu2[2] * S22) / (
            S22 + wavelength_variance_local
        )
        x0 = matrix.col((mu1[0], mu1[1])) + S12 * (1 / S22) * (z0 - mu1[2])
        mu3 = matrix.col((x0[0], x0[1], z0))

        # Rotate the sigma and mean vector back
        Sigma3 = R.transpose() * Sigma3 * R
        s3 = R.transpose() * mu3
        q3 = s3 - s0

        # Compute the mean wavelength and diffracted beam vector
        wavelength3 = -2 * s0.normalize().dot(q3) / (q3.dot(q3))
        mu_s1 = s0.normalize() / wavelength3 + q3

        # Compute vector to approximate spread in diffracted vectors
        v = -mu_s1 / (s0.dot(q3))

        # Compute the vector to approximate the spread in diffracted vectors
        # v = -0.5*(2*q3/(s0.dot(q3)) - q3.dot(q3)/s0.dot(q3)**2*s0)

        # Compute the mean diffracted vector
        # mu_s1 = (v.dot(q3))*s0 + q3

        # try:
        #   assert (mu_s1 - mu_s1_2).length() < 1e-10
        #   assert (v - v_2).length() < 1e-10
        # except Exception:
        #   print tuple(s0)
        #   print tuple(q3)
        #   print s0.dot(q3)
        #   print tuple(mu_s1), tuple(mu_s1_2)
        #   print tuple(v), tuple(v_2)
        #   raise

        # Compute the sigma of s1
        sigma_q = Sigma3
        sigma_vs = (v.transpose() * sigma_q * v)[0] * s0.dot(s0)
        eigen_values_v = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, sigma_vs))
        sigma_v = R0.transpose() * eigen_values_v * R0
        sigma_vq = s0 * (
            v.transpose() * (sigma_q + q3 * q3.transpose())
            - (v.dot(q3)) * q3.transpose()
        )
        sigma_qv = sigma_vq.transpose()
        sigma_s1 = sigma_q + sigma_v + sigma_vq + sigma_qv

        # Rotate into the coordinate system along mean s1
        # x_axis = mu_s1.cross(s0).normalize()
        # z_axis = mu_s1.normalize()
        # y_axis = z_axis.cross(x_axis).normalize()
        # R2 = matrix.sqr(
        #   x_axis.elems +
        #   y_axis.elems +
        #   z_axis.elems)

        # Rotate the sigma and marginalize along that direction
        Sigma4 = R * sigma_s1 * R.transpose()
        Sigma5 = matrix.sqr((Sigma4[0], Sigma4[1], Sigma4[3], Sigma4[4]))
        # Sigma5_inv = Sigma5.inverse()
        mu4 = R * mu_s1
        mu5 = matrix.col((mu4[0], mu4[1]))

        Scale = matrix.sqr((1, 0, 0, 1)) * (s0.length() / mu_s1.length())

        Sigma6 = Scale.transpose() * Sigma5 * Scale
        mu6 = Scale * mu5

        # # Marginal sigma and mean
        # Sigma_marginal = matrix.sqr((
        #   Sigma3[0], Sigma3[1],
        #   Sigma3[3], Sigma3[4]))
        # mu_marginal = matrix.col((mu3[0], mu3[1]))

        # Compute rotated mu
        # mu = R*s2
        # assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
        # mu1 = matrix.col((mu[0], mu[1]))
        # mu2 = mu[2]

        # Compute conditional mean and covariance
        # mubar = matrix.col((0,0)) + S12*(1/S22)*(s0.length()-s2.length())
        # Sbar = S11 - S12*(1/S22)*S21
        # S2 = S22

        Sbar = Sigma6
        mubar = mu6
        S2 = S22 + wavelength_variance_local
        # print tuple(Sbar), tuple(Sobs)
        # Compute the log likelihood
        try:
            Sobs = matrix.sqr(Sobs)
            B1 = ctot * log(S2)
            B2 = ctot * (1 / S2) * (s0.length() - s2.length()) ** 2
            A1 = log(Sbar.determinant()) * ctot
            A2 = ctot * (Sbar.inverse() * Sobs).trace()
            A3 = (
                ctot
                * (
                    Sbar.inverse() * ((xbar - mubar) * (xbar - mubar).transpose())
                ).trace()
            )
            lnL += -0.5 * (A1 + A2 + A3 + (B1 + B2))
        except Exception:
            raise
            lnL += -1e15

    if wavelength_variance > 0.10 ** 2:
        lnL += -1e15
    # if max(spot_covariance[0], spot_covariance[4], spot_covariance[8]) > 1e-5:
    #   lnL += -1e15

    # v = 1000000
    # p = 3
    # P = matrix.sqr((
    #   1e-7, 0, 0,
    #   0, 1e-7, 0,
    #   0, 0, 1e-7))
    # X = spot_covariance
    # l = -(v+p+1)*log(X.determinant())/2 -  0.5*(P*X.inverse()).trace()
    # print tuple(spot_covariance), lnL, l
    # lnL += l

    # v = 1
    # wv0 = 1e-7
    # if (wavelength_variance > 1e-20):
    #   lnL += -(v*wv0)/(2*wavelength_variance) - (1+v/2)*log(wavelength_variance)
    # else:
    #   lnL += -(v*wv0)/(2*1e-20) - (1+v/2)*log(1e-20)

    # print wavelength_variance, ("%.2e " * 9) % tuple(spot_covariance), lnL
    # print tuple(sigma), lnL
    return lnL


class Target(object):
    def __init__(self, s0, s2_list, xbar_list, ctot_list, Sobs_list):
        self.s0 = s0
        self.s2_list = s2_list
        self.xbar_list = xbar_list
        self.ctot_list = ctot_list
        self.Sobs_list = Sobs_list

    def __call__(self, params, *args):
        params = flex.double(params)
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


def plot_function(func):

    X = []
    Y = []
    params = (
        -0.0018922555580294159,
        -0.003960342417917075,
        0.000802512703902016,
        -0.011212605823672782,
        0.0015804959560581232,
        0.0011843711796837466,
        8.205550283434665e-08,
    )
    # params = (0.9912606131407175,
    #           0.8721521410328605, -0.0005427636101740637,
    #           -0.4755941521897193, 6.978595033167225e-05, -0.0005654927049233458,
    #           0.3446544460324457)
    # params = (0.0004,
    #           1e-7, -0.0005427636101740637,
    #           1e-5, 6.978595033167225e-05, -0.0005654927049233458,
    #           0.3446544460324457)
    # wavelength_param = tan(sqrt(0.000111616168007) *pi / (2*0.05))

    minx = -0.00001
    maxx = 0.00001
    for i in range(200):
        param = minx + i * (maxx - minx) / 200.0
        values = flex.double(params)
        values[0] = param
        # values = flex.double((
        #   paramssqrt(1e-7),
        #   0, sqrt(2e-7),
        #   0, 0, sqrt(param),
        #   wavelength_param ))

        x = param  # (0.05*atan(param) / (pi / 2))**2
        y = func(values)
        print(i, x, y)

        X.append(x)
        Y.append(y)

    from matplotlib import pylab

    pylab.plot(X, Y)
    pylab.show()
    exit(0)


def tst_ideal():

    numpy.random.seed(100)

    # The beam vector
    s0 = matrix.col((0, 0, 1.03))

    A = matrix.sqr(
        (
            0.009186272656211013,
            -0.0009923913997518274,
            0.009833912026596807,
            0.004663921486320244,
            -0.0005161976250487236,
            -0.019341482302357254,
            0.007239709464264245,
            0.012542003882687962,
            -1.793473438470342e-05,
        )
    )

    q0 = A * matrix.col((0, 0, 0))
    q1 = A * matrix.col((1, 0, 0))
    q2 = A * matrix.col((0, 1, 0))
    q3 = A * matrix.col((0, 0, 1))
    d1 = (q0 - q1).length()
    d2 = (q0 - q2).length()
    d3 = (q0 - q3).length()
    print(min((d1 / 3.0) ** 2, (d2 / 3.0) ** 2, (d3 / 3.0) ** 2))

    origin = matrix.col((-217.87240000000003, 227.8312, 305.0))
    fast_axis = matrix.col((1, 0, 0))
    slow_axis = matrix.col((0, -1, 0))

    D = (
        matrix.sqr(fast_axis.elems + slow_axis.elems + origin.elems)
        .transpose()
        .inverse()
    )

    # The covariance matrix
    spot_covariance = matrix.sqr((1e-7, 0, 0, 0, 2e-7, 0, 0, 0, 3e-7))

    wavelength_variance = 0.05 ** 2  # 1e-4 # variance

    # The number of reflections
    N = 100

    # Generate a load of reflections
    s2_list, ctot_list, xbar_list, Sobs_list = generate_with_wavelength_spread2(
        D, s0, spot_covariance, wavelength_variance, N=N
    )

    wavelength_param = sqrt(1e-3)  # 0#sqrt(0.05**2)#tan(sqrt(1e-12) *pi / (2*0.05))

    # plot_function(Target(s0, s2_list, xbar_list, ctot_list, Sobs_list))

    # Starting values for simplex
    # values = flex.double((
    #   sqrt(1e-7),
    #   0, sqrt(1e-7),
    #   0, 0, sqrt(1e-7)))
    # values = flex.double((
    #   sqrt(1e-7),sqrt(1e-7),sqrt(1e-7), 0, 0, 0, wavelength_param))
    values = flex.double(
        (sqrt(1e-7), 0, sqrt(1e-7), 0, 0, sqrt(1e-7), wavelength_param)
    )
    # values = flex.double([wavelength_param])
    # sqrt(0.01**2)))
    # offset = flex.double(
    #   [sqrt(1e-7)  for v in values])

    target = Target(s0, s2_list, xbar_list, ctot_list, Sobs_list)

    def callback(x):
        spot_covariance = matrix_from_params(x)
        wavelength_variance = x[6] ** 2  # (0.05*atan(x[6]) / (pi / 2))**2
        lnL = -target(x)
        print(wavelength_variance, ("%.2e " * 9) % tuple(spot_covariance), lnL)

    # Do the simplex optimization
    result = minimize(
        target,
        values,
        #    method="Nelder-Mead")
        method="BFGS",
        callback=callback,
    )
    params = flex.double(result.x)
    # optimizer = SimpleSimplex(
    #   values,
    #   offset,
    #   Target(
    #     s0,
    #     s2_list,
    #     xbar_list,
    #     ctot_list,
    #     Sobs_list),
    #   2000)
    # params = optimizer.get_solution()

    # Create the covariance matrix
    spot_covariance = matrix_from_params(params)

    print(spot_covariance)

    wavelength_variance = params[6] ** 2  # (0.05*atan(params[6]) / (pi / 2))**2
    print(wavelength_variance)

    # expected = matrix.sqr((
    #   9.91047018199e-07, -1.98078253593e-09, 2.27093231797e-09,
    #   -1.98078253593e-09, 1.98335548957e-06, 1.88051940862e-08,
    #   2.27093231797e-09, 1.88051940862e-08, 2.99885951955e-06))
    # assert all(1e6*abs(a-b) < 1e-7 for a, b in zip(sigma, expected))
    print(tuple(params))

    print("OK")


if __name__ == "__main__":
    tst_ideal()
