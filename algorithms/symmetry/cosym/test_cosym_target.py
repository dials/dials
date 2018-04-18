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

  for weights in [None, 'count', 'standard_error']:

    t = target.Target(datasets,
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
    assert g.all_approx_equal_relatively(g_fd, relative_error=1e-4)

    c = t.curvatures(x)
    c_fd = t.curvatures_fd(x, eps=1e-3)
    assert c.all_approx_equal_relatively(c_fd, relative_error=0.5e-1)

    M = engine.lbfgs_with_curvs(
      target=t, coords=x,
      verbose=False,
    )
    t.compute_functional(x)
    # check functional has decreased and gradients are approximately zero
    f, g = t.compute_functional_and_gradients(x)
    g_fd = t.compute_gradients_fd(x)
    assert f < f0
    assert g.all_approx_equal(0, 1e-4)
    assert g_fd.all_approx_equal(0, 1e-4)

