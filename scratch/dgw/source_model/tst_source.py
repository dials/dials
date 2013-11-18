from __future__ import division
from dials.scratch.dgw.source_model import source
from scitbx import matrix

if __name__ == '__main__':

  # Really dumb test

  import random
  from libtbx.test_utils import approx_equal

  s0 = source(matrix.col.random(3, 0.5, 1.5))
  curr_lambda = s0.get_wavelength()
  s0.set_s0(matrix.col.random(3, 0.5, 1.5))
  assert approx_equal(s0.get_wavelength(), curr_lambda)

  print "OK"
