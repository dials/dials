from __future__ import annotations

import numpy as np
import pytest

from dials.algorithms.scaling.tof.lambda_norm import Observations, fit_lambda_norm
from dials.array_family import flex
from dials_tof_scaling_ext import (
    GPRWavelengthScale,
    LambdaNormVarianceModel,
    leave_one_out_samples,
)


def spectrum(lam):
    return 3.0 * np.exp(-(((lam - 2.4) / 0.9) ** 2)) + 0.3


def synthetic(rng, w=spectrum, n_groups=4000, multiplicity=6, factors=(1.0,)):
    F = rng.gamma(1.0, 2000.0, n_groups)
    lam = rng.uniform(1.5, 4.0, (n_groups, multiplicity))
    expt = rng.integers(0, len(factors), (n_groups, multiplicity))
    mean = F[:, None] * w(lam) * np.asarray(factors)[expt]
    var = mean + 200.0
    I = rng.normal(mean, np.sqrt(var))
    offsets = np.arange(0, n_groups * multiplicity + 1, multiplicity)
    return Observations(
        intensity=I.ravel(),
        variance=var.ravel(),
        wavelength=lam.ravel(),
        experiment=expt.ravel(),
        group_offsets=offsets,
    )


@pytest.fixture
def params():
    from dials.command_line.lambda_norm import phil_scope

    return phil_scope.extract()


def test_variance_model_recovers_counting_term():
    rng = np.random.default_rng(1)
    lam = rng.uniform(1, 4, 5000)
    I = rng.gamma(1.0, 1000.0, 5000)
    vm = LambdaNormVarianceModel()
    vm.fit(flex.double(lam), flex.double(I), flex.double(0.5 * I + 50.0))
    assert list(vm.b_values()) == pytest.approx([0.5] * 8, rel=1e-6)
    assert vm.variance(2.0, 100.0, 100.0, 300.0) == pytest.approx(200.0)

    vm.fit(flex.double(lam), flex.double(I), flex.double(5000, 50.0))
    assert list(vm.b_values()) == pytest.approx([0.0] * 8, abs=1e-12)


def test_variance_model_too_few_observations():
    vm = LambdaNormVarianceModel()
    vm.fit(flex.double(10, 2.0), flex.double(10, 1.0), flex.double(10, 1.0))
    assert vm.b(2.0) == 0
    assert vm.variance(2.0, 5.0, 3.0, 50.0) == 3.0


def test_leave_one_out_matches_numpy():
    rng = np.random.default_rng(2)
    I = rng.uniform(10, 100, 7)
    var = rng.uniform(5, 20, 7)
    lam = rng.uniform(1, 3, 7)
    corr = rng.uniform(0.5, 2, 7)
    offsets = [0, 1, 4, 7]  # the singleton group gives no samples
    s = leave_one_out_samples(
        flex.double(I),
        flex.double(var),
        flex.double(lam),
        flex.double(corr),
        flex.size_t(offsets),
        LambdaNormVarianceModel(),
        False,
    )
    assert list(s.index) == [1, 2, 3, 4, 5, 6]
    for k, i in enumerate(s.index):
        g = 1 if i < 4 else 2
        mates = [j for j in range(offsets[g], offsets[g + 1]) if j != i]
        w = corr[mates] ** 2 / var[mates]
        M = np.sum(w * I[mates] / corr[mates]) / np.sum(w)
        assert s.rho[k] == pytest.approx(I[i] / M)
        assert s.u[k] == pytest.approx(M / var[i])


def test_fit_recovers_spectrum(params):
    obs = synthetic(np.random.default_rng(0))
    params.reference_wavelength = 2.4
    gpr = fit_lambda_norm(obs, params).gpr
    lam = np.linspace(1.6, 3.9, 12)
    w = np.array(gpr.scale(flex.double(lam)))
    assert w == pytest.approx(spectrum(lam) / spectrum(2.4), rel=0.02)
    assert gpr.scale(flex.double([2.4]))[0] == pytest.approx(1.0, rel=1e-3)


def test_fit_with_experiment_factors(params):
    obs = synthetic(np.random.default_rng(3), factors=(0.5, 1.0, 2.0))
    params.reference_wavelength = 2.4
    result = fit_lambda_norm(obs, params)
    k = result.experiment_factors
    assert k / k[1] == pytest.approx([0.5, 1.0, 2.0], rel=0.02)
    lam = np.linspace(1.6, 3.9, 12)
    w = np.array(result.gpr.scale(flex.double(lam)))
    assert w == pytest.approx(spectrum(lam) / spectrum(2.4), rel=0.02)


def test_fit_flat_spectrum(params):
    obs = synthetic(np.random.default_rng(4), w=lambda lam: np.ones_like(lam))
    gpr = fit_lambda_norm(obs, params).gpr
    w = np.array(gpr.scale(flex.double(np.linspace(1.5, 4.0, 20))))
    assert w == pytest.approx(1.0, abs=0.02)


def test_held_at_range_ends(params):
    obs = synthetic(np.random.default_rng(5))
    gpr = fit_lambda_norm(obs, params).gpr
    w = gpr.scale(flex.double([0.5, gpr.lambda_min(), gpr.lambda_max(), 10.0]))
    assert w[0] == w[1]
    assert w[2] == w[3]


def test_too_few_samples():
    gpr = GPRWavelengthScale(50, 0.0, False, 0.0, 0.0)
    assert not gpr.fit(
        flex.double(np.linspace(1, 3, 10)),
        flex.double(10, 1.0),
        flex.double(10, 1.0),
        2.0,
    )
    assert "too few" in gpr.message()
    assert list(gpr.scale(flex.double([1.0, 2.0]))) == [1.0, 1.0]
