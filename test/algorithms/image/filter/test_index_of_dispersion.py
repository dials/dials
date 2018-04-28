from __future__ import absolute_import, division, print_function

import random

import pytest

def test():
  from dials.algorithms.image.filter import index_of_dispersion_filter
  from scitbx.array_family import flex

  # Create an image
  image = flex.random_double(2000 * 2000)
  image.reshape(flex.grid(2000, 2000))
  mask = flex.random_bool(2000 * 2000, 0.99).as_int()
  mask.reshape(flex.grid(2000, 2000))

  # Calculate the summed area table
  mask2 = mask.deep_copy()
  index_of_dispersion_filter = index_of_dispersion_filter(image, mask2, (3, 3), 2)
  mean = index_of_dispersion_filter.mean()
  var = index_of_dispersion_filter.sample_variance()
  index_of_dispersion = index_of_dispersion_filter.index_of_dispersion()

  # For a selection of random points, ensure that the value is the
  # sum of the area under the kernel
  eps = 1e-7
  for i in range(10000):
    i = random.randint(10, 1990)
    j = random.randint(10, 1990)
    m1 = mean[j,i]
    v1 = var[j,i]
    f1 = index_of_dispersion[j,i]
    p = image[j-3:j+4,i-3:i+4]
    m = mask[j-3:j+4,i-3:i+4]
    if mask[j,i] == 0:
      m2 = 0.0
      v2 = 0.0
      f2 = 1.0
    else:
      p = flex.select(p, flags=m)
      mv = flex.mean_and_variance(flex.double(p))
      m2 = mv.mean()
      v2 = mv.unweighted_sample_variance()
      f2 = v2 / m2
    assert m1 == pytest.approx(m2, abs=eps)
    assert v1 == pytest.approx(v2, abs=eps)
    assert f1 == pytest.approx(f2, abs=eps)
