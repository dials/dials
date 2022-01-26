from __future__ import annotations

import numpy as np
import pytest

from cctbx import sgtbx
from scitbx.array_family import flex

from dials.algorithms.symmetry.cosym import engine, target
from dials.algorithms.symmetry.cosym._generate_test_data import generate_test_data


@pytest.mark.parametrize("space_group", ["P2", "P3", "P6", "R3:h", "I23"])
def test_cosym_target(space_group):
    datasets, expected_reindexing_ops = generate_test_data(
        space_group=sgtbx.space_group_info(symbol=space_group).group(), sample_size=50
    )

    intensities = datasets[0]
    dataset_ids = np.zeros(intensities.size() * len(datasets))
    for i, d in enumerate(datasets[1:]):
        i += 1
        intensities = intensities.concatenate(d, assert_is_similar_symmetry=False)
        dataset_ids[i * d.size() : (i + 1) * d.size()] = np.full(d.size(), i, dtype=int)

    for weights in [None, "count", "standard_error"]:
        print(weights)
        t = target.Target(intensities, dataset_ids, weights=weights)
        m = len(t.sym_ops)
        n = len(datasets)
        assert t.dim == m
        assert t.rij_matrix.shape == (n * m, n * m)
        # x = np.random.rand(n * m * t.dim)
        x = flex.random_double(n * m * t.dim).as_numpy_array()
        f0 = t.compute_functional(x)
        g = t.compute_gradients(x)
        g_fd = t.compute_gradients_fd(x)
        np.testing.assert_allclose(g, g_fd, rtol=2e-3)
        c = t.curvatures(x)
        c_fd = t.curvatures_fd(x, eps=1e-3)
        assert list(c) == pytest.approx(c_fd, rel=0.8e-1)

        if weights == "count":
            # Absolute upper limit on weights
            assert t.wij_matrix.max() <= datasets[0].size()

        minimizer = engine.lbfgs_with_curvs(target=t, coords=x)
        # check functional has decreased and gradients are approximately zero
        f = t.compute_functional(minimizer.coords)
        g = t.compute_gradients(minimizer.coords)
        g_fd = t.compute_gradients_fd(minimizer.coords)
        assert f < f0
        assert pytest.approx(g, abs=1e-3) == [0] * len(g)
        assert pytest.approx(g_fd, abs=1e-3) == [0] * len(g)
