from __future__ import absolute_import, division, print_function

import pytest

from cctbx import sgtbx
from scitbx.array_family import flex

from dials.algorithms.symmetry.cosym.generate_test_data import generate_test_data
from dials.algorithms.symmetry.cosym import engine
from dials.algorithms.symmetry.cosym import target


@pytest.mark.parametrize('space_group', ['P2', 'P3', 'P6', 'R3:h', 'I23'])
def test_cosym_target(space_group):
  datasets, expected_reindexing_ops = generate_test_data(
    space_group=sgtbx.space_group_info(symbol=space_group).group())

  intensities = datasets[0]
  dataset_ids = flex.double(intensities.size(), 0)
  for i, d in enumerate(datasets[1:]):
    intensities = intensities.concatenate(
      d, assert_is_similar_symmetry=False)
    dataset_ids.extend(flex.double(d.size(), i+1))

  for weights in [None, 'count', 'standard_error']:

    t = target.Target(intensities,
                      dataset_ids,
                      weights=weights,
                      )
    m = len(t.get_sym_ops())
    n = len(datasets)
    assert t.dim == m
    assert t.rij_matrix.all() == (n*m, n*m)
    x = flex.random_double(n * m * t.dim)
    x_orig = x.deep_copy()
    f0, g = t.compute_functional_and_gradients(x)
    g_fd = t.compute_gradients_fd(x)
    assert list(g) == pytest.approx(g_fd, rel=1e-3)

    c = t.curvatures(x)
    c_fd = t.curvatures_fd(x, eps=1e-3)
    assert list(c) == pytest.approx(c_fd, rel=0.5e-1)

    M = engine.lbfgs_with_curvs(
      target=t, coords=x,
      verbose=False,
    )
    t.compute_functional(x)
    # check functional has decreased and gradients are approximately zero
    f, g = t.compute_functional_and_gradients(x)
    g_fd = t.compute_gradients_fd(x)
    assert f < f0
    assert pytest.approx(g, abs=1e-3) == [0]*len(g)
    assert pytest.approx(g_fd, abs=1e-3) == [0]*len(g)
