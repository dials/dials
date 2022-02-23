from __future__ import annotations

from random import uniform

import pytest

from scitbx import matrix


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
        return (-(1 / 12) * (A + E) + (4 / 3) * (B + D) - (5 / 2) * C) / h**2
    else:
        A = func(x - h, y - h)
        B = func(x - h, y)
        C = func(x, y - h)
        D = func(x, y)
        E = func(x, y + h)
        F = func(x + h, y)
        G = func(x + h, y + h)
        return (A - B - C + 2 * D - E - F + G) / (2 * h**2)


def generate_data():

    from random import seed

    seed(0)

    s2 = matrix.col(
        (uniform(0, 1), uniform(0, 1), uniform(0, 1))
    ).normalize() * uniform(0.9, 1.1)

    s0 = matrix.col((0, 0, 1))

    b1, b2, b3, b4, b5, b6 = (
        uniform(0.5, 3.5),
        uniform(0.0, 0.5),
        uniform(0.5, 3.5),
        uniform(0.0, 0.5),
        uniform(0.0, 0.5),
        uniform(0.5, 3.5),
    )

    M = matrix.sqr((b1, 0, 0, b2, b3, 0, b4, b5, b6))

    sigma = M * M.transpose()

    # ctot = randint(100, 1000)

    return sigma, s0, s2


def generate_testdata():
    for i in range(10):
        sigma, s0, (b1, b2, b3) = generate_data()
        yield sigma, s0, b1, b2, b3


def ds2_db(b1, b2, b3):
    return [matrix.col((1, 0, 0)), matrix.col((0, 1, 0)), matrix.col((0, 0, 1))]


def compute_s2(b1, b2, b3):
    return matrix.col((b1, b2, b3))


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_epsilon(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_epsilon(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        return s0.length() - s2.length()

    def f1(x):
        return compute_epsilon(x, b2, b3)

    def f2(x):
        return compute_epsilon(b1, x, b3)

    def f3(x):
        return compute_epsilon(b1, b2, x)

    h = 0.001

    depdb1_num = first_derivative(f1, b1, h)
    depdb2_num = first_derivative(f2, b2, h)
    depdb3_num = first_derivative(f3, b3, h)

    def compute_dep(s2, ds2):
        return -(1 / s2.length()) * s2.dot(ds2)

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    depdb1_cal = compute_dep(s2, ds2[0])
    depdb2_cal = compute_dep(s2, ds2[1])
    depdb3_cal = compute_dep(s2, ds2[2])

    assert abs(depdb1_num - depdb1_cal) < 1e-7
    assert abs(depdb2_num - depdb2_cal) < 1e-7
    assert abs(depdb3_num - depdb3_cal) < 1e-7


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_mubar(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_mubar(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        epsilon = s0.length() - s2.length()
        return S12 * (1 / S22) * epsilon

    def f1(x):
        return compute_mubar(x, b2, b3)

    def f2(x):
        return compute_mubar(b1, x, b3)

    def f3(x):
        return compute_mubar(b1, b2, x)

    h = 0.001

    dmubar_db1_num = first_derivative(f1, b1, h)
    dmubar_db2_num = first_derivative(f2, b2, h)
    dmubar_db3_num = first_derivative(f3, b3, h)

    def compute_depsilon(s2, ds2):
        return -(1 / s2.length()) * s2.dot(ds2)

    def compute_dmubar(s2, ds2):
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        return S12 * (1 / S22) * compute_depsilon(s2, ds2)

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    dmubar_db1_cal = compute_dmubar(s2, ds2[0])
    dmubar_db2_cal = compute_dmubar(s2, ds2[1])
    dmubar_db3_cal = compute_dmubar(s2, ds2[2])

    assert abs(dmubar_db1_num - dmubar_db1_cal) < 1e-7
    assert abs(dmubar_db2_num - dmubar_db2_cal) < 1e-7
    assert abs(dmubar_db3_num - dmubar_db3_cal) < 1e-7


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_e1(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_e1(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        return s2.cross(s0).normalize()

    def f1(x):
        return compute_e1(x, b2, b3)

    def f2(x):
        return compute_e1(b1, x, b3)

    def f3(x):
        return compute_e1(b1, b2, x)

    h = 0.001

    de1db1_num = first_derivative(f1, b1, h)
    de1db2_num = first_derivative(f2, b2, h)
    de1db3_num = first_derivative(f3, b3, h)

    def compute_de1(s2, ds2):
        e1 = s2.cross(s0)
        d1 = ds2.cross(s0)
        A = d1 / e1.length()
        B = e1 * d1.dot(e1) / (e1.length() ** 3)
        return A - B

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    de1db1_cal = compute_de1(s2, ds2[0])
    de1db2_cal = compute_de1(s2, ds2[1])
    de1db3_cal = compute_de1(s2, ds2[2])

    assert abs(de1db1_num - de1db1_cal) < 1e-7
    assert abs(de1db2_num - de1db2_cal) < 1e-7
    assert abs(de1db3_num - de1db3_cal) < 1e-7


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_e2(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_e2(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        e1 = s2.cross(s0).normalize()
        return s2.cross(e1).normalize()

    def f1(x):
        return compute_e2(x, b2, b3)

    def f2(x):
        return compute_e2(b1, x, b3)

    def f3(x):
        return compute_e2(b1, b2, x)

    h = 0.001

    de2db1_num = first_derivative(f1, b1, h)
    de2db2_num = first_derivative(f2, b2, h)
    de2db3_num = first_derivative(f3, b3, h)

    def compute_de2(s2, ds2):
        e1 = s2.cross(s0)
        e2 = s2.cross(e1)
        d2 = ds2 * s2.dot(s0) + s2 * ds2.dot(s0) - 2 * s0 * s2.dot(ds2)
        A = d2 / e2.length()
        B = e2 * d2.dot(e2) / (e2.length() ** 3)
        return A - B

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    de2db1_cal = compute_de2(s2, ds2[0])
    de2db2_cal = compute_de2(s2, ds2[1])
    de2db3_cal = compute_de2(s2, ds2[2])

    assert abs(de2db1_num - de2db1_cal) < 1e-7
    assert abs(de2db2_num - de2db2_cal) < 1e-7
    assert abs(de2db3_num - de2db3_cal) < 1e-7


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_e3(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_e3(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        return s2.normalize()

    def f1(x):
        return compute_e3(x, b2, b3)

    def f2(x):
        return compute_e3(b1, x, b3)

    def f3(x):
        return compute_e3(b1, b2, x)

    h = 0.001

    de3db1_num = first_derivative(f1, b1, h)
    de3db2_num = first_derivative(f2, b2, h)
    de3db3_num = first_derivative(f3, b3, h)

    def compute_de3(s2, ds2):
        e3 = s2
        d3 = ds2
        A = d3 / e3.length()
        B = e3 * d3.dot(e3) / (e3.length() ** 3)
        return A - B

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    de3db1_cal = compute_de3(s2, ds2[0])
    de3db2_cal = compute_de3(s2, ds2[1])
    de3db3_cal = compute_de3(s2, ds2[2])

    assert abs(de3db1_num - de3db1_cal) < 1e-7
    assert abs(de3db2_num - de3db2_cal) < 1e-7
    assert abs(de3db3_num - de3db3_cal) < 1e-7


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_s1(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_s1(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        epsilon = s0.length() - s2.length()
        mubar = S12 * (1 / S22) * epsilon

        e1 = s2.cross(s0).normalize()
        e2 = s2.cross(e1).normalize()
        e3 = s2.normalize()
        R = matrix.sqr(e1.elems + e2.elems + e3.elems)
        s1 = R.transpose() * matrix.col((mubar[0], mubar[1], s0.length()))
        return s1

    def f1(x):
        return compute_s1(x, b2, b3)

    def f2(x):
        return compute_s1(b1, x, b3)

    def f3(x):
        return compute_s1(b1, b2, x)

    h = 0.001

    db1_num = first_derivative(f1, b1, h)
    db2_num = first_derivative(f2, b2, h)
    db3_num = first_derivative(f3, b3, h)

    def compute_depsilon(s2, ds2):
        return -(1 / s2.length()) * s2.dot(ds2)

    def compute_dmubar(s2, ds2):
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        return S12 * (1 / S22) * compute_depsilon(s2, ds2)

    def compute_R(s2):
        e1 = s2.cross(s0).normalize()
        e2 = s2.cross(e1).normalize()
        e3 = s2.normalize()
        R = matrix.sqr(e1.elems + e2.elems + e3.elems)
        return R

    def compute_de1(s2, ds2):
        e1 = s2.cross(s0)
        d1 = ds2.cross(s0)
        A = d1 / e1.length()
        B = e1 * d1.dot(e1) / (e1.length() ** 3)
        return A - B

    def compute_de2(s2, ds2):
        e1 = s2.cross(s0)
        e2 = s2.cross(e1)
        d2 = ds2 * s2.dot(s0) + s2 * ds2.dot(s0) - 2 * s0 * s2.dot(ds2)
        A = d2 / e2.length()
        B = e2 * d2.dot(e2) / (e2.length() ** 3)
        return A - B

    def compute_de3(s2, ds2):
        e3 = s2
        d3 = ds2
        A = d3 / e3.length()
        B = e3 * d3.dot(e3) / (e3.length() ** 3)
        return A - B

    def compute_dR(s2, ds2):
        de1 = compute_de1(s2, ds2)
        de2 = compute_de2(s2, ds2)
        de3 = compute_de3(s2, ds2)
        return matrix.sqr(de1.elems + de2.elems + de3.elems)

    def compute_ds1(s2, ds2):
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        R = compute_R(s2)
        dR = compute_dR(s2, ds2)
        epsilon = s0.length() - s2.length()
        mubar = S12 * (1 / S22) * epsilon
        dmubar = compute_dmubar(s2, ds2)
        v = matrix.col((mubar[0], mubar[1], s0.length()))
        A = -R.transpose() * dR * R.transpose() * v
        B = R.transpose() * matrix.col((dmubar[0], dmubar[1], 0))
        return A + B

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    db1_cal = compute_ds1(s2, ds2[0])
    db2_cal = compute_ds1(s2, ds2[1])
    db3_cal = compute_ds1(s2, ds2[2])

    assert abs(db1_num - db1_cal) < 1e-7
    assert abs(db2_num - db2_cal) < 1e-7
    assert abs(db3_num - db3_cal) < 1e-7


@pytest.mark.parametrize("sigma,s0,b1,b2,b3", generate_testdata())
def test_derivative_of_f(sigma, s0, b1, b2, b3):

    ds2 = ds2_db(b1, b2, b3)

    def compute_f(b1, b2, b3):
        s2 = compute_s2(b1, b2, b3)
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        epsilon = s0.length() - s2.length()
        mubar = S12 * (1 / S22) * epsilon

        e1 = s2.cross(s0).normalize()
        e2 = s2.cross(e1).normalize()
        e3 = s2.normalize()
        R = matrix.sqr(e1.elems + e2.elems + e3.elems)
        s1 = R.transpose() * matrix.col((mubar[0], mubar[1], s0.length()))

        D = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

        v = D * s1
        assert v[2] > 0
        X = v[0] / v[2]
        Y = v[1] / v[2]

        xcal = matrix.col((X, Y))
        xobs = matrix.col((0, 0))
        return (xobs - xcal).transpose() * (xobs - xcal)

    def f1(x):
        return compute_f(x, b2, b3)

    def f2(x):
        return compute_f(b1, x, b3)

    def f3(x):
        return compute_f(b1, b2, x)

    h = 0.001

    db1_num = first_derivative(f1, b1, h)
    db2_num = first_derivative(f2, b2, h)
    db3_num = first_derivative(f3, b3, h)

    def compute_depsilon(s2, ds2):
        return -(1 / s2.length()) * s2.dot(ds2)

    def compute_dmubar(s2, ds2):
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        return S12 * (1 / S22) * compute_depsilon(s2, ds2)

    def compute_R(s2):
        e1 = s2.cross(s0).normalize()
        e2 = s2.cross(e1).normalize()
        e3 = s2.normalize()
        R = matrix.sqr(e1.elems + e2.elems + e3.elems)
        return R

    def compute_de1(s2, ds2):
        e1 = s2.cross(s0)
        d1 = ds2.cross(s0)
        A = d1 / e1.length()
        B = e1 * d1.dot(e1) / (e1.length() ** 3)
        return A - B

    def compute_de2(s2, ds2):
        e1 = s2.cross(s0)
        e2 = s2.cross(e1)
        d2 = ds2 * s2.dot(s0) + s2 * ds2.dot(s0) - 2 * s0 * s2.dot(ds2)
        A = d2 / e2.length()
        B = e2 * d2.dot(e2) / (e2.length() ** 3)
        return A - B

    def compute_de3(s2, ds2):
        e3 = s2
        d3 = ds2
        A = d3 / e3.length()
        B = e3 * d3.dot(e3) / (e3.length() ** 3)
        return A - B

    def compute_dR(s2, ds2):
        de1 = compute_de1(s2, ds2)
        de2 = compute_de2(s2, ds2)
        de3 = compute_de3(s2, ds2)
        return matrix.sqr(de1.elems + de2.elems + de3.elems)

    def compute_ds1(s2, ds2):
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        R = compute_R(s2)
        dR = compute_dR(s2, ds2)
        epsilon = s0.length() - s2.length()
        mubar = S12 * (1 / S22) * epsilon
        dmubar = compute_dmubar(s2, ds2)
        v = matrix.col((mubar[0], mubar[1], s0.length()))
        A = -R.transpose() * dR * R.transpose() * v
        B = R.transpose() * matrix.col((dmubar[0], dmubar[1], 0))
        return A + B

    def compute_s1(s2):
        S12 = matrix.col((sigma[2], sigma[5]))
        S22 = sigma[8]
        epsilon = s0.length() - s2.length()
        mubar = S12 * (1 / S22) * epsilon

        e1 = s2.cross(s0).normalize()
        e2 = s2.cross(e1).normalize()
        e3 = s2.normalize()
        R = matrix.sqr(e1.elems + e2.elems + e3.elems)
        s1 = R.transpose() * matrix.col((mubar[0], mubar[1], s0.length()))
        return s1

    def compute_df(s2, ds2):
        s1 = compute_s1(s2)
        ds1 = compute_ds1(s2, ds2)
        D = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
        v = D * s1
        # dv = D * ds1
        assert v[2] > 0
        X = v[0] / v[2]
        Y = v[1] / v[2]
        xcal = matrix.col((X, Y))
        dx = (v[2] * ds1[0] - v[0] * ds1[2]) / v[2] ** 2
        dy = (v[2] * ds1[1] - v[1] * ds1[2]) / v[2] ** 2
        xobs = matrix.col((0, 0))
        dxy = matrix.col((dx, dy))
        return -2 * (xobs - xcal).transpose() * dxy

    s2 = compute_s2(b1, b2, b3)
    ds2 = ds2_db(b1, b2, b3)

    db1_cal = compute_df(s2, ds2[0])
    db2_cal = compute_df(s2, ds2[1])
    db3_cal = compute_df(s2, ds2[2])

    assert abs(db1_num - db1_cal) < 1e-7
    assert abs(db2_num - db2_cal) < 1e-7
    assert abs(db3_num - db3_cal) < 1e-7
