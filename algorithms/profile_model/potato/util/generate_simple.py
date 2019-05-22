from math import exp, sqrt

from numpy.random import multivariate_normal, normal, uniform

from scitbx import matrix

from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)
from dials.array_family import flex


def generate_simple(s0, sigma, N=100):
    """
    Generate a list of normally distributed observations

    """
    s2_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    # Loop through the list
    for k in range(N):

        # Compute position in reciprocal space of the centre of the rlp
        s2_direction = matrix.col(
            (uniform(0, 1), uniform(0, 1), uniform(0, 1))
        ).normalize()

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, s2_direction)
        sigmap = R * sigma * R.transpose()
        s2_magnitude = normal(s0.length(), sqrt(sigmap[8]))
        s2 = s2_direction * s2_magnitude

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Partition the matrix
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]

        # Compute the conditional distribution
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]
        z = s0.length()
        sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21
        mu_bar = mu1 + sigma12 * (z - mu2) / sigma22

        # Compute the scale factor and intensity
        # scale = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # Simulate some observations
        points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))

        # Compute the observed mean for each observation
        ctot = 0
        xbar = matrix.col((0, 0))
        for x in points:
            xbar += matrix.col(x)

        ctot = len(points)

        xbar /= ctot

        # Compute the observed covariance for each observation
        Sobs = matrix.sqr((0, 0, 0, 0))
        for x in points:
            x = matrix.col(x)
            Sobs += (x - xbar) * (x - xbar).transpose()
        Sobs /= ctot

        s2_list.append(s2)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return s2_list, ctot_list, xbar_list, Sobs_list


def generate_simple_binned(s0, sigma, N=100):
    """
    Generate a list of normally distributed observations

    """
    s2_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    # Loop through the list
    for k in range(N):

        # Compute position in reciprocal space of the centre of the rlp
        s2_direction = matrix.col(
            (uniform(0, 1), uniform(0, 1), uniform(0, 1))
        ).normalize()

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, s2_direction)
        sigmap = R * sigma * R.transpose()
        s2_magnitude = normal(s0.length(), sqrt(sigmap[8]))
        s2 = s2_direction * s2_magnitude

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Partition the matrix
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]

        # Compute the conditional distribution
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]
        z = s0.length()
        sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21
        mu_bar = mu1 + sigma12 * (z - mu2) / sigma22

        # Compute the scale factor and intensity
        # scale = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # Simulate some observations
        a = 7.5
        b = 15.0 / (12 * sqrt(sigma_bar[0]))
        c = 15.0 / (12 * sqrt(sigma_bar[3]))
        D = flex.double(flex.grid(15, 15))
        points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))
        for x, y in points:
            i = int(a + b * x)
            j = int(a + c * y)
            if j >= 0 and j < D.all()[0] and i >= 0 and i < D.all()[1]:
                D[j, i] += 1

        # Compute the observed mean for each observation
        ctot = 0
        xbar = matrix.col((0, 0))
        for j in range(D.all()[0]):
            for i in range(D.all()[1]):
                ctot += D[j, i]
                x = matrix.col(((i + 0.5 - a) / b, (j + 0.5 - a) / c))
                xbar += x * D[j, i]

        if ctot <= 0:
            continue

        xbar /= ctot

        # Compute the observed covariance for each observation
        Sobs = matrix.sqr((0, 0, 0, 0))
        for j in range(D.all()[0]):
            for i in range(D.all()[1]):
                x = matrix.col(((i + 0.5 - a) / b, (j + 0.5 - a) / c))
                Sobs += (x - xbar) * (x - xbar).transpose() * D[j, i]
        Sobs /= ctot

        s2_list.append(s2)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return s2_list, ctot_list, xbar_list, Sobs_list


def generate_from_reflections(s0, sigma, reflections):
    """
    Generate a list of normally distributed observations

    """
    s2_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    # Loop through the list
    for k in range(len(reflections)):

        # Compute position in reciprocal space of the centre of the rlp
        s2 = matrix.col(reflections["s2"][k])

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, s2)
        sigmap = R * sigma * R.transpose()

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Partition the matrix
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]

        # Compute the conditional distribution
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]
        z = s0.length()
        sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21
        mu_bar = mu1 + sigma12 * (z - mu2) / sigma22

        # Perform rejection sampling to get a normally distributed set of
        # reflections
        P = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        R = uniform(0, 1.0)
        if P < R:
            continue

        # Compute the scale factor and intensity
        # scale = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # Simulate some observations
        points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))

        # Compute the observed mean for each observation
        ctot = 0
        xbar = matrix.col((0, 0))
        for x in points:
            xbar += matrix.col(x)

        ctot = len(points)

        xbar /= ctot

        # Compute the observed covariance for each observation
        Sobs = matrix.sqr((0, 0, 0, 0))
        for x in points:
            x = matrix.col(x)
            Sobs += (x - xbar) * (x - xbar).transpose()
        Sobs /= ctot

        s2_list.append(s2)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return s2_list, ctot_list, xbar_list, Sobs_list


def generate_from_reflections2(A, s0, sigma, reflections):
    """
    Generate a list of normally distributed observations

    """
    sp_list = []
    h_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    # Loop through the list
    for k in range(len(reflections)):

        # Compute position in reciprocal space of the centre of the rlp
        h = matrix.col(reflections["miller_index"][k])
        r = matrix.sqr(A) * h
        s2 = s0 + r
        sp = s2

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, sp)
        sigmap = R * sigma * R.transpose()

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Partition the matrix
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]

        # Compute the conditional distribution
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]
        z = s0.length()
        sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21
        mu_bar = mu1 + sigma12 * (z - mu2) / sigma22

        # Perform rejection sampling to get a normally distributed set of
        # reflections
        P = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        R = uniform(0, 1.0)
        if P < R:
            continue

        # Compute the scale factor and intensity
        # scale = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # Simulate some observations
        points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))

        # Compute the observed mean for each observation
        ctot = 0
        xbar = matrix.col((0, 0))
        for x in points:
            xbar += matrix.col(x)

        ctot = len(points)

        xbar /= ctot

        # Compute the observed covariance for each observation
        Sobs = matrix.sqr((0, 0, 0, 0))
        for x in points:
            x = matrix.col(x)
            Sobs += (x - xbar) * (x - xbar).transpose()
        Sobs /= ctot

        sp_list.append(sp)
        h_list.append(h)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return sp_list, h_list, ctot_list, xbar_list, Sobs_list


def generate_from_reflections_binned(s0, sigma, reflections):
    """
    Generate a list of normally distributed observations

    """
    s2_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    # Loop through the list
    for k in range(len(reflections)):

        # Compute position in reciprocal space of the centre of the rlp
        s2 = matrix.col(reflections["s2"][k])

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, s2)
        sigmap = R * sigma * R.transpose()

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Partition the matrix
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]

        # Compute the conditional distribution
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]
        z = s0.length()
        sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21
        mu_bar = mu1 + sigma12 * (z - mu2) / sigma22

        # Perform rejection sampling to get a normally distributed set of
        # reflections
        P = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        R = uniform(0, 1.0)
        if P < R:
            continue

        # Compute the scale factor and intensity
        # scale = exp(-0.5 * (z - mu2) ** 2 / sigma22)
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # Simulate some observations
        a = 7.5
        b = 15.0 / (12 * sqrt(sigma_bar[0]))
        c = 15.0 / (12 * sqrt(sigma_bar[3]))
        D = flex.double(flex.grid(15, 15))
        points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))
        for x, y in points:
            i = int(a + b * x)
            j = int(a + c * y)
            if j >= 0 and j < D.all()[0] and i >= 0 and i < D.all()[1]:
                D[j, i] += 1

        # Compute the observed mean for each observation
        ctot = 0
        xbar = matrix.col((0, 0))
        for j in range(D.all()[0]):
            for i in range(D.all()[1]):
                ctot += D[j, i]
                x = matrix.col(((i + 0.5 - a) / b, (j + 0.5 - a) / c))
                xbar += x * D[j, i]

        if ctot <= 0:
            continue

        xbar /= ctot

        # Compute the observed covariance for each observation
        Sobs = matrix.sqr((0, 0, 0, 0))
        for j in range(D.all()[0]):
            for i in range(D.all()[1]):
                x = matrix.col(((i + 0.5 - a) / b, (j + 0.5 - a) / c))
                Sobs += (x - xbar) * (x - xbar).transpose() * D[j, i]
        Sobs /= ctot

        s2_list.append(s2)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return s2_list, ctot_list, xbar_list, Sobs_list


def generate_with_wavelength_spread(s0, sigma_spot, sigma_wavelength, N=100):
    """
    Generate a list of normally distributed observations

    """
    s2_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    # Loop through the list
    for k in range(N):

        # Compute position in reciprocal space of the centre of the rlp
        s2_direction = matrix.col(
            (uniform(0, 1), uniform(0, 1), uniform(0, 1))
        ).normalize()

        # Compute the local spread in energy
        q0 = s2_direction.normalize() * s0.length() - s0
        sigma_wavelength_local = sqrt(sigma_wavelength) * q0.dot(q0) / 2

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, s2_direction)
        sigmap = R * sigma_spot * R.transpose()
        s2_magnitude = normal(
            s0.length(), sqrt(sigmap[8] + sigma_wavelength_local ** 2)
        )
        # s2_magnitude = normal(s0.length(), sqrt(sigmap[8]))
        s2 = s2_direction * s2_magnitude

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Apply the wavelength spread to the sigma
        Sigma1_inv = sigmap.inverse()
        Sigma2_inv = matrix.sqr(
            (0, 0, 0, 0, 0, 0, 0, 0, 1 / sigma_wavelength_local ** 2)
        )
        Sigma3_inv = Sigma1_inv + Sigma2_inv
        Sigma3 = Sigma3_inv.inverse()

        # Apply the wavelength spread to the mean
        mu1 = mu
        mu2 = matrix.col((0, 0, s0.length()))
        mu3 = Sigma3 * (Sigma1_inv * mu1 + Sigma2_inv * mu2)

        # Marginal sigma and mean
        Sigma_marginal = matrix.sqr((Sigma3[0], Sigma3[1], Sigma3[3], Sigma3[4]))
        mu_marginal = matrix.col((mu3[0], mu3[1]))

        # Compute the scale factor and intensity
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # # Simulate some observations
        # points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))

        # # Compute the observed mean for each observation
        # ctot = 0
        # xbar = matrix.col((0, 0))
        # for x in points:
        #   xbar += matrix.col(x)

        # ctot = len(points)

        # xbar /= ctot

        # # Compute the observed covariance for each observation
        # Sobs = matrix.sqr((0, 0, 0, 0))
        # for x in points:
        #   x = matrix.col(x)
        #   Sobs += (x-xbar)*(x-xbar).transpose()
        # Sobs /= ctot

        ctot = I
        xbar = mu_marginal
        Sobs = Sigma_marginal

        s2_list.append(s2)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return s2_list, ctot_list, xbar_list, Sobs_list


def generate_with_wavelength_spread2(
    D, s0, spot_covariance, wavelength_variance, N=100
):
    """
    Generate a list of normally distributed observations

    """
    s2_list = []
    ctot_list = []
    xbar_list = []
    Sobs_list = []

    D_inv = D.inverse()

    # Loop through the list
    for k in range(N):

        # Compute position in reciprocal space of the centre of the rlp
        s2_direction = matrix.col(
            (uniform(0, 1), uniform(0, 1), uniform(1, 1))
        ).normalize()

        # Compute the local spread in energy
        q0 = s2_direction.normalize() * s0.length() - s0
        wavelength_variance_local = wavelength_variance * (q0.dot(q0) / 2) ** 2

        # Rotate the covariance matrix
        R = compute_change_of_basis_operation(s0, s2_direction)
        sigmap = R * spot_covariance * R.transpose()
        sigma_spread = sigmap[8] + wavelength_variance_local
        s2_magnitude = normal(s0.length(), sqrt(sigma_spread))
        # s2_magnitude = normal(s0.length(), sqrt(sigmap[8]))
        s2 = s2_direction * s2_magnitude

        # resolution = 1.0 / (2.0 * s0.length() * sin(0.5 * s0.angle(s2)))

        # print resolution

        # Rotate to get mu
        mu = R * s2
        assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

        # Partition matrix
        S11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        S12 = matrix.col((sigmap[2], sigmap[5]))
        S21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        S22 = sigmap[8]

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

        # Apply the wavelength spread to the mean
        mu1 = mu
        mu2 = matrix.col((0, 0, s0.length()))
        z0 = (mu1[2] * wavelength_variance_local + mu2[2] * S22) / (
            S22 + wavelength_variance_local
        )
        x0 = matrix.col((mu1[0], mu1[1])) + S12 * (1 / S22) * (z0 - mu1[2])
        mu3 = matrix.col((x0[0], x0[1], z0))

        # Apply the wavelength spread to the sigma
        # Sigma1_inv = sigmap.inverse()
        # Sigma2_inv = matrix.sqr((
        #   0, 0, 0,
        #   0, 0, 0,
        #   0, 0, 1/wavelength_variance_local))
        # Sigma3_inv = Sigma1_inv + Sigma2_inv
        # Sigma3 = Sigma3_inv.inverse()

        # # Apply the wavelength spread to the mean
        # mu1 = mu
        # mu2 = matrix.col((0, 0, s0.length()))
        # mu3 = Sigma3*(Sigma1_inv*mu1 + Sigma2_inv*mu2)

        # Rotate the sigma and mean vector back
        Sigma3 = R.transpose() * Sigma3 * R
        s3 = R.transpose() * mu3
        q3 = s3 - s0

        # Compute the scale factor and intensity
        I = uniform(50, 1000)
        if I <= 1:
            continue

        # Simulate some observations
        points = multivariate_normal(q3, Sigma3.as_list_of_lists(), int(I))

        X = []
        Y = []

        # Generate points
        for p in points:

            p = matrix.col(p)

            # Compute wavelength and diffracting vector for point
            wavelength = -2.0 * s0.normalize().dot(p) / (p.dot(p))
            s0w = s0.normalize() / wavelength
            s1w = s0w + p
            assert abs(s1w.length() - s0w.length()) < 1e-10

            # Do the ray projection onto the detector
            v = D * s1w
            assert v[2] > 0
            x = v[0] / v[2]
            y = v[1] / v[2]
            X.append(x)
            Y.append(y)

        # Compute observed mean on detector
        # xmean = sum(X) / len(X)
        # ymean = sum(Y) / len(Y)

        # Map onto coordinate system and compute mean
        ctot = len(X)
        xbar = matrix.col((0, 0))
        # cs = CoordinateSystem2d(s0, s2)
        for x, y in zip(X, Y):
            s = (D_inv * matrix.col((x, y, 1))).normalize() * s0.length()
            e = R * (s - s2.normalize() * s0.length())
            e1, e2 = e[0], e[1]
            # e1, e2 = cs.from_beam_vector(s)
            xbar += matrix.col((e1, e2))
        xbar /= ctot

        # Compute variance
        Sobs = matrix.sqr((0, 0, 0, 0))
        for x, y in zip(X, Y):
            s = (D_inv * matrix.col((x, y, 1))).normalize() * s0.length()
            e = R * (s - s2.normalize() * s0.length())
            e1, e2 = e[0], e[1]
            # e1, e2 = cs.from_beam_vector(s)
            x = matrix.col((e1, e2))
            Sobs += (x - xbar) * (x - xbar).transpose()
        Sobs /= ctot

        # print tuple(xbar), tuple(Sobs)

        s2_list.append(s2)
        ctot_list.append(ctot)
        xbar_list.append(xbar)
        Sobs_list.append(list(Sobs))

    return s2_list, ctot_list, xbar_list, Sobs_list
