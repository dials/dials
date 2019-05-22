from __future__ import division, print_function

from math import log
from random import randint, uniform

import pytest

from scitbx import matrix

from dials.algorithms.profile_model.potato.parameterisation import (
    Simple6MosaicityParameterisation,
)
from dials.algorithms.profile_model.potato.refiner import RefinerData


def first_derivative(func, x, h):
    return (-func(x + 2 * h) + 8 * func(x + h) - 8 * func(x - h) + func(x - 2 * h)) / (
        12 * h
    )


def second_derivative(func, x, y=None, h=None):
    if y is None:
        A = func(x + 2 * h)
        B = func(x + h)
        C = func(x)
        D = func(x - h)
        E = func(x - 2 * h)
        return (-(1 / 12) * (A + E) + (4 / 3) * (B + D) - (5 / 2) * C) / h ** 2
    else:
        A = func(x - h, y - h)
        B = func(x - h, y)
        C = func(x, y - h)
        D = func(x, y)
        E = func(x, y + h)
        F = func(x + h, y)
        G = func(x + h, y + h)
        return (A - B - C + 2 * D - E - F + G) / (2 * h ** 2)


def generate_data():

    from random import seed

    seed(0)

    s2 = matrix.col(
        (uniform(0, 1), uniform(0, 1), uniform(0, 1))
    ).normalize() * uniform(0.9, 1.1)

    s0 = matrix.col((0, 0, 1))

    T = matrix.sqr((uniform(0.5, 2.5), 0, uniform(0, 0.1), uniform(0.5, 2.5)))
    S = T * T.transpose()

    X = matrix.col((0.01, -0.03))

    # b1, b2, b3, b4, b5, b6 = 1, 0.1, 2, 0.2, 0.3, 3
    b1, b2, b3, b4, b5, b6 = (
        uniform(0.5, 3.5),
        uniform(0.0, 0.5),
        uniform(0.5, 3.5),
        uniform(0.0, 0.5),
        uniform(0.0, 0.5),
        uniform(0.5, 3.5),
    )

    ctot = randint(100, 1000)

    return (b1, b2, b3, b4, b5, b6), s0, s2, ctot, X, S


def generate_testdata():

    for i in range(10):

        (b1, b2, b3, b4, b5, b6), s0, s2, ctot, mobs, Sobs = generate_data()

        parameterisation = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
        reflection_model = RefinerData(parameterisation, s0, s2, ctot, mobs, Sobs)

        yield reflection_model


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_dSdb_22(reflection_model):

    (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
    R = reflection_model.R
    # mu2 = reflection_model.s2.length()
    # r = reflection_model.s0.length()
    # Sobs = reflection_model.sobs
    # ctot = reflection_model.ctot

    def compute_sigma22(b1, b2, b3, b4, b5, b6):
        model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
        sigma = model.sigma()
        sigmap = R * sigma * R.transpose()
        sigma22 = sigmap[8]
        return sigma22

    def f1(x):
        return compute_sigma22(x, b2, b3, b4, b5, b6)

    def f2(x):
        return compute_sigma22(b1, x, b3, b4, b5, b6)

    def f3(x):
        return compute_sigma22(b1, b2, x, b4, b5, b6)

    def f4(x):
        return compute_sigma22(b1, b2, b3, x, b5, b6)

    def f5(x):
        return compute_sigma22(b1, b2, b3, b4, x, b6)

    def f6(x):
        return compute_sigma22(b1, b2, b3, b4, b5, x)

    h = 0.001

    dSdb1_22_num = first_derivative(f1, b1, h)
    dSdb2_22_num = first_derivative(f2, b2, h)
    dSdb3_22_num = first_derivative(f3, b3, h)
    dSdb4_22_num = first_derivative(f4, b4, h)
    dSdb5_22_num = first_derivative(f5, b5, h)
    dSdb6_22_num = first_derivative(f6, b6, h)

    dSdb1_22_cal = reflection_model.marginal.first_derivatives()[0]
    dSdb2_22_cal = reflection_model.marginal.first_derivatives()[1]
    dSdb3_22_cal = reflection_model.marginal.first_derivatives()[2]
    dSdb4_22_cal = reflection_model.marginal.first_derivatives()[3]
    dSdb5_22_cal = reflection_model.marginal.first_derivatives()[4]
    dSdb6_22_cal = reflection_model.marginal.first_derivatives()[5]

    assert abs(dSdb1_22_num - dSdb1_22_cal) < 1e-7
    assert abs(dSdb2_22_num - dSdb2_22_cal) < 1e-7
    assert abs(dSdb3_22_num - dSdb3_22_cal) < 1e-7
    assert abs(dSdb4_22_num - dSdb4_22_cal) < 1e-7
    assert abs(dSdb5_22_num - dSdb5_22_cal) < 1e-7
    assert abs(dSdb6_22_num - dSdb6_22_cal) < 1e-7

    # print 'OK'


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_dm_bar_db(reflection_model):

    (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
    R = reflection_model.R
    s0 = reflection_model.s0
    mu2 = reflection_model.s2.length()
    # r = reflection_model.s0.length()
    # Sobs = reflection_model.sobs
    # ctot = reflection_model.ctot
    epsilon = s0.length() - mu2

    def compute_mu_bar(b1, b2, b3, b4, b5, b6):
        model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
        sigma = model.sigma()
        sigmap = R * sigma * R.transpose()
        # sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        # sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]
        return sigma12 * (1 / sigma22) * epsilon

    def f1(x):
        return compute_mu_bar(x, b2, b3, b4, b5, b6)

    def f2(x):
        return compute_mu_bar(b1, x, b3, b4, b5, b6)

    def f3(x):
        return compute_mu_bar(b1, b2, x, b4, b5, b6)

    def f4(x):
        return compute_mu_bar(b1, b2, b3, x, b5, b6)

    def f5(x):
        return compute_mu_bar(b1, b2, b3, b4, x, b6)

    def f6(x):
        return compute_mu_bar(b1, b2, b3, b4, b5, x)

    h = 0.001

    dmubar_db1_num = first_derivative(f1, b1, h)
    dmubar_db2_num = first_derivative(f2, b2, h)
    dmubar_db3_num = first_derivative(f3, b3, h)
    dmubar_db4_num = first_derivative(f4, b4, h)
    dmubar_db5_num = first_derivative(f5, b5, h)
    dmubar_db6_num = first_derivative(f6, b6, h)

    dmubar_db1_cal = reflection_model.conditional.first_derivatives_of_mean()[0]
    dmubar_db2_cal = reflection_model.conditional.first_derivatives_of_mean()[1]
    dmubar_db3_cal = reflection_model.conditional.first_derivatives_of_mean()[2]
    dmubar_db4_cal = reflection_model.conditional.first_derivatives_of_mean()[3]
    dmubar_db5_cal = reflection_model.conditional.first_derivatives_of_mean()[4]
    dmubar_db6_cal = reflection_model.conditional.first_derivatives_of_mean()[5]

    assert all(abs(a - b) < 1e-7 for a, b in zip(dmubar_db1_num, dmubar_db1_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dmubar_db2_num, dmubar_db2_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dmubar_db3_num, dmubar_db3_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dmubar_db4_num, dmubar_db4_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dmubar_db5_num, dmubar_db5_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dmubar_db6_num, dmubar_db6_cal))


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_dS_bar_db(reflection_model):

    (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
    R = reflection_model.R
    # mu2 = reflection_model.s2.length()
    # r = reflection_model.s0.length()
    # Sobs = reflection_model.sobs
    # ctot = reflection_model.ctot

    def compute_sigma_bar(b1, b2, b3, b4, b5, b6):
        model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
        sigma = model.sigma()
        sigmap = R * sigma * R.transpose()
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]
        return sigma11 - sigma12 * (1 / sigma22) * sigma21

    def f1(x):
        return compute_sigma_bar(x, b2, b3, b4, b5, b6)

    def f2(x):
        return compute_sigma_bar(b1, x, b3, b4, b5, b6)

    def f3(x):
        return compute_sigma_bar(b1, b2, x, b4, b5, b6)

    def f4(x):
        return compute_sigma_bar(b1, b2, b3, x, b5, b6)

    def f5(x):
        return compute_sigma_bar(b1, b2, b3, b4, x, b6)

    def f6(x):
        return compute_sigma_bar(b1, b2, b3, b4, b5, x)

    h = 0.001

    dS_bar_db1_num = first_derivative(f1, b1, h)
    dS_bar_db2_num = first_derivative(f2, b2, h)
    dS_bar_db3_num = first_derivative(f3, b3, h)
    dS_bar_db4_num = first_derivative(f4, b4, h)
    dS_bar_db5_num = first_derivative(f5, b5, h)
    dS_bar_db6_num = first_derivative(f6, b6, h)

    dS_bar_db1_cal = reflection_model.conditional.first_derivatives_of_sigma()[0]
    dS_bar_db2_cal = reflection_model.conditional.first_derivatives_of_sigma()[1]
    dS_bar_db3_cal = reflection_model.conditional.first_derivatives_of_sigma()[2]
    dS_bar_db4_cal = reflection_model.conditional.first_derivatives_of_sigma()[3]
    dS_bar_db5_cal = reflection_model.conditional.first_derivatives_of_sigma()[4]
    dS_bar_db6_cal = reflection_model.conditional.first_derivatives_of_sigma()[5]

    assert all(abs(a - b) < 1e-7 for a, b in zip(dS_bar_db1_num, dS_bar_db1_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dS_bar_db2_num, dS_bar_db2_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dS_bar_db3_num, dS_bar_db3_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dS_bar_db4_num, dS_bar_db4_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dS_bar_db5_num, dS_bar_db5_cal))
    assert all(abs(a - b) < 1e-7 for a, b in zip(dS_bar_db6_num, dS_bar_db6_cal))

    # print 'OK'


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_dLdb(reflection_model):

    (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
    R = reflection_model.R
    s0 = reflection_model.s0
    mu2 = reflection_model.s2.length()
    r = reflection_model.s0.length()
    mobs = reflection_model.mobs
    Sobs = reflection_model.sobs
    ctot = reflection_model.ctot

    def compute_L(b1, b2, b3, b4, b5, b6):
        model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
        sigma = model.sigma()
        sigmap = R * sigma * R.transpose()
        sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
        sigma12 = matrix.col((sigmap[2], sigmap[5]))
        sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
        sigma22 = sigmap[8]

        z = s0.length()
        mubar = sigma12 * (1 / sigma22) * (z - mu2)
        sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21

        d = r - mu2
        c_d = mubar - mobs
        A = log(sigma22)
        B = (1 / sigma22) * d ** 2
        C = log(sigma_bar.determinant()) * ctot
        D = (sigma_bar.inverse() * ctot * Sobs).trace()
        E = (sigma_bar.inverse() * ctot * c_d * c_d.transpose()).trace()
        return -0.5 * (A + B + C + D + E)

    def f1(x):
        return compute_L(x, b2, b3, b4, b5, b6)

    def f2(x):
        return compute_L(b1, x, b3, b4, b5, b6)

    def f3(x):
        return compute_L(b1, b2, x, b4, b5, b6)

    def f4(x):
        return compute_L(b1, b2, b3, x, b5, b6)

    def f5(x):
        return compute_L(b1, b2, b3, b4, x, b6)

    def f6(x):
        return compute_L(b1, b2, b3, b4, b5, x)

    h = 0.001

    dLdb1_num = first_derivative(f1, b1, h)
    dLdb2_num = first_derivative(f2, b2, h)
    dLdb3_num = first_derivative(f3, b3, h)
    dLdb4_num = first_derivative(f4, b4, h)
    dLdb5_num = first_derivative(f5, b5, h)
    dLdb6_num = first_derivative(f6, b6, h)

    dL = reflection_model.first_derivatives()

    dLdb1_cal = dL[0]
    dLdb2_cal = dL[1]
    dLdb3_cal = dL[2]
    dLdb4_cal = dL[3]
    dLdb5_cal = dL[4]
    dLdb6_cal = dL[5]

    assert abs(dLdb1_num - dLdb1_cal) < 1e-7
    assert abs(dLdb2_num - dLdb2_cal) < 1e-7
    assert abs(dLdb3_num - dLdb3_cal) < 1e-7
    assert abs(dLdb4_num - dLdb4_cal) < 1e-7
    assert abs(dLdb5_num - dLdb5_cal) < 1e-7
    assert abs(dLdb6_num - dLdb6_cal) < 1e-7

    # print 'OK'


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_d2S_dbij(reflection_model):
    def test_for_index(i, j):
        (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
        R = reflection_model.R
        # mu2 = reflection_model.s2.length()
        # r = reflection_model.s0.length()
        # Sobs = reflection_model.sobs
        # ctot = reflection_model.ctot

        def compute_S(b1, b2, b3, b4, b5, b6):
            model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
            sigma = model.sigma()
            sigmap = R * sigma * R.transpose()
            # sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
            # sigma12 = matrix.col((sigmap[2], sigmap[5]))
            # sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
            sigma22 = sigmap[8]
            return sigma22

        def f1(x):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            return compute_S(*params)

        def f2(x, y):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            params[j] = y
            return compute_S(*params)

        h = 0.001

        x = (b1, b2, b3, b4, b5, b6)[i]
        y = (b1, b2, b3, b4, b5, b6)[j]

        if i == j:
            d2Sdb_22_num = second_derivative(f1, x=x, h=h)
        else:
            d2Sdb_22_num = second_derivative(f2, x=x, y=y, h=h)

        d2Sdb_22_cal = reflection_model.marginal.second_derivatives()[i, j]

        assert abs(d2Sdb_22_num - d2Sdb_22_cal) < 1e-7

    for j in range(6):
        for i in range(6):
            test_for_index(i, j)


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_d2S_bar_dbij(reflection_model):
    def test_for_index(i, j):
        (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
        R = reflection_model.R
        # mu2 = reflection_model.s2.length()
        # r = reflection_model.s0.length()
        # Sobs = reflection_model.sobs
        # ctot = reflection_model.ctot

        def compute_Sbar(b1, b2, b3, b4, b5, b6):
            model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
            sigma = model.sigma()
            sigmap = R * sigma * R.transpose()
            sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
            sigma12 = matrix.col((sigmap[2], sigmap[5]))
            sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
            sigma22 = sigmap[8]
            return sigma11 - sigma12 * (1 / sigma22) * sigma21

        def f1(x):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            return compute_Sbar(*params)

        def f2(x, y):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            params[j] = y
            return compute_Sbar(*params)

        h = 0.001

        x = (b1, b2, b3, b4, b5, b6)[i]
        y = (b1, b2, b3, b4, b5, b6)[j]

        if i == j:
            d2S_bar_db_num = second_derivative(f1, x=x, h=h)
        else:
            d2S_bar_db_num = second_derivative(f2, x=x, y=y, h=h)

        d2S_bar_db_cal = reflection_model.conditional.second_derivatives_of_sigma()[i][
            j
        ]

        assert all(abs(a - b) < 1e-5 for a, b in zip(d2S_bar_db_num, d2S_bar_db_cal))

    for j in range(6):
        for i in range(6):
            test_for_index(i, j)


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_d2m_bar_dbij(reflection_model):
    def test_for_index(i, j):
        (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()

        R = reflection_model.R
        s0 = reflection_model.s0
        mu2 = reflection_model.s2.length()
        # r = reflection_model.s0.length()
        # Sobs = reflection_model.sobs
        # ctot = reflection_model.ctot
        epsilon = s0.length() - mu2

        def compute_mu_bar(b1, b2, b3, b4, b5, b6):
            model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
            sigma = model.sigma()
            sigmap = R * sigma * R.transpose()
            # sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
            sigma12 = matrix.col((sigmap[2], sigmap[5]))
            # sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
            sigma22 = sigmap[8]
            return sigma12 * (1 / sigma22) * epsilon

        def f1(x):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            return compute_mu_bar(*params)

        def f2(x, y):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            params[j] = y
            return compute_mu_bar(*params)

        h = 0.001

        x = (b1, b2, b3, b4, b5, b6)[i]
        y = (b1, b2, b3, b4, b5, b6)[j]

        if i == j:
            d2m_bar_db_num = second_derivative(f1, x=x, h=h)
        else:
            d2m_bar_db_num = second_derivative(f2, x=x, y=y, h=h)

        d2m_bar_db_cal = reflection_model.conditional.second_derivatives_of_mean()[i][j]

        assert all(abs(a - b) < 1e-5 for a, b in zip(d2m_bar_db_num, d2m_bar_db_cal))

    for j in range(6):
        for i in range(6):
            test_for_index(i, j)


@pytest.mark.parametrize("reflection_model", generate_testdata())
def test_d2L_dbij(reflection_model):
    def test_for_index(i, j):
        (b1, b2, b3, b4, b5, b6) = reflection_model.model.parameters()
        R = reflection_model.R
        mu2 = reflection_model.s2.length()
        r = reflection_model.s0.length()
        mobs = reflection_model.mobs
        Sobs = reflection_model.sobs
        ctot = reflection_model.ctot
        s0 = reflection_model.s0

        def compute_L(b1, b2, b3, b4, b5, b6):
            model = Simple6MosaicityParameterisation((b1, b2, b3, b4, b5, b6))
            sigma = model.sigma()
            sigmap = R * sigma * R.transpose()
            sigma11 = matrix.sqr((sigmap[0], sigmap[1], sigmap[3], sigmap[4]))
            sigma12 = matrix.col((sigmap[2], sigmap[5]))
            sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
            sigma22 = sigmap[8]

            z = s0.length()
            mubar = sigma12 * (1 / sigma22) * (z - mu2)
            sigma_bar = sigma11 - sigma12 * (1 / sigma22) * sigma21

            d = r - mu2
            c_d = mubar - mobs
            A = log(sigma22)
            B = (1 / sigma22) * d ** 2
            C = log(sigma_bar.determinant()) * ctot
            D = (sigma_bar.inverse() * ctot * Sobs).trace()
            E = (sigma_bar.inverse() * ctot * c_d * c_d.transpose()).trace()
            return -0.5 * (A + B + C + D + E)

        def f1(x):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            return compute_L(*params)

        def f2(x, y):
            params = [b1, b2, b3, b4, b5, b6]
            params[i] = x
            params[j] = y
            return compute_L(*params)

        h = 0.001

        x = (b1, b2, b3, b4, b5, b6)[i]
        y = (b1, b2, b3, b4, b5, b6)[j]

        if i == j:
            d2Ldb_num = second_derivative(f1, x=x, h=h)
        else:
            d2Ldb_num = second_derivative(f2, x=x, y=y, h=h)

        d2Ldb_cal = reflection_model.second_derivatives()[i, j]

        # print d2Ldb_cal, d2Ldb_num

        assert abs(d2Ldb_num - d2Ldb_cal) < 1e-3

    for j in range(6):
        for i in range(6):
            test_for_index(i, j)
