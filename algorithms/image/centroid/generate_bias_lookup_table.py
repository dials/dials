"""
Code to generate lookup values for algorithms/image/centroid/bias.h.
"""
from __future__ import absolute_import, division, print_function


def sum_of_erf(mu, sigma, N=1000):
    """
    Compute the sum of erf term

    :param mu: The Gaussian mean
    :param sigma: The Gaussian sigma
    :param N: The number of iterations in the sum

    """
    from math import sqrt, erf

    sum1 = 0
    sum2 = N * erf((N + 1 - mu) / (sqrt(2) * sigma))
    sum3 = N * erf((-N - mu) / (sqrt(2) * sigma))
    for i in range(-N, N):
        sum1 += erf((i + 1 - mu) / (sqrt(2) * sigma))
    return -0.5 * (sum1 + sum2 + sum3)


def compute_normal_bias_sq(sigma, N1=1000, N2=1000):
    """
    Compute the bias for a Gaussian profile

    :param sigma: The Gaussian sigma
    :param N1: The number of iterations for the sum_of_erf
    :param N2: The number of divisions for the integral

    """

    # If sigma is zero then we know the bias is 1/12.0
    if sigma == 0:
        return 1.0 / 12.0

    # The function to integrate
    def function(mu):
        return (sum_of_erf(mu, sigma, N1) - mu + 0.5) ** 2

    # Perform the numerical integration
    fa = function(0.0)
    fb = function(1.0)
    sum_f = 0.0
    for i in range(1, N2):
        sum_f += function(i * 1.0 / N2)
    return (1.0 / N2) * (fa / 2.0 + fb / 2.0 + sum_f)


def compute_lookup_table(max_sigma=0.5, N1=1000, N2=1000, N3=50):
    """
    Compute a lookup table of bias for Gaussian sigma values

    :param max_sigma: The maximum sigma to compute
    :param N1: The number of iterations for the sum_of_erf
    :param N2: The number of division for the integral
    :param N3: The number of elements in the table

    """
    sigma = []
    bias_sq = []
    for i in range(N3):
        s = i * max_sigma / N3
        b = compute_normal_bias_sq(s, N1, N2)
        sigma.append(s)
        bias_sq.append(b)

    return sigma, bias_sq


if __name__ == "__main__":

    sigma, bias_sq = compute_lookup_table()

    for s, b in zip(sigma, bias_sq):
        print("%0.2f %0.7f" % (s, b))
