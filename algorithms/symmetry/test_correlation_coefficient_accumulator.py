from __future__ import absolute_import, division, print_function

import pytest

from scitbx.array_family import flex
from dials.algorithms.symmetry.determine_space_group import CorrelationCoefficientAccumulator

def test_correlation_coefficient_accumulator():
  x = flex.double((1,2,3))
  cc_accumulator = CorrelationCoefficientAccumulator(x, x)
  corr = flex.linear_correlation(x, x)
  assert cc_accumulator.coefficient() == pytest.approx(corr.coefficient())
  assert cc_accumulator.n() == x.size()
  assert cc_accumulator.coefficient() == pytest.approx(1.0)

  # compare with flex.linear_correlation()
  n = 100
  x = flex.random_double(n)
  y = flex.random_double(n)
  cc_accumulator = CorrelationCoefficientAccumulator(x, y)
  corr = flex.linear_correlation(x, y)
  assert cc_accumulator.coefficient() == pytest.approx(corr.coefficient())
  assert cc_accumulator.n() == x.size()

  # test CorrelationCoefficientAccumulator.accumulate()
  x_all = x
  y_all = y
  x_ = flex.random_double(n)
  y_ = flex.random_double(n)
  cc_accumulator.accumulate(x_, y_)
  x_all.extend(x_)
  y_all.extend(y_)
  corr = flex.linear_correlation(x_all, y_all)
  assert cc_accumulator.coefficient() == pytest.approx(corr.coefficient())
  assert cc_accumulator.n() == x_all.size()

  # test CorrelationCoefficientAccumulator += other
  x_ = flex.random_double(n)
  y_ = flex.random_double(n)
  cc_accumulator += CorrelationCoefficientAccumulator(x_, y_)
  x_all.extend(x_)
  y_all.extend(y_)
  corr = flex.linear_correlation(x_all, y_all)
  assert cc_accumulator.coefficient() == pytest.approx(corr.coefficient())
  assert cc_accumulator.n() == x_all.size()

