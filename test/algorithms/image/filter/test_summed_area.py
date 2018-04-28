from __future__ import absolute_import, division, print_function

import random

import pytest

def test():
  from dials.algorithms.image.filter import summed_area
  from scitbx.array_family import flex

  # Create an image
  image = flex.random_double(2000 * 2000)
  image.reshape(flex.grid(2000, 2000))

  # Calculate the summed area table
  sa = summed_area(image, (3, 3))

  # For a selection of random points, ensure that the value is the
  # sum of the area under the kernel
  for i in range(10000):
    i = random.randint(10, 1990)
    j = random.randint(10, 1990)
    v = sa[j,i]
    e = flex.sum(image[j-3:j+4,i-3:i+4])
    assert e == pytest.approx(v, abs=1e-7)
