#!/usr/bin/env python

from __future__ import absolute_import, division
from math import pi
import random

from libtbx.test_utils import approx_equal
from scitbx import matrix

from dxtbx.model import BeamFactory
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisation
from dials.algorithms.refinement.refinement_helpers \
    import get_fd_gradients, random_param_shift

if __name__ == '__main__':

  # make a random beam vector and parameterise it
  bf = BeamFactory()
  s0 = bf.make_beam(matrix.col.random(3, 0.5, 1.5), wavelength=1.2)
  s0p = BeamParameterisation(s0)

  # Let's do some basic tests. First, can we change parameter values and
  # update the modelled vector s0?
  s0_old = matrix.col(s0.get_s0())
  s0p.set_param_vals([1000*0.1, 1000*0.1, 0.8])
  assert(approx_equal(matrix.col(s0.get_s0()).angle(s0_old), 0.1413033))
  assert(approx_equal(matrix.col(s0.get_s0()).length(), 0.8))

  # random initial orientations and wavelengths with a random parameter shifts
  attempts = 1000
  failures = 0
  for i in range(attempts):

    # make a random beam vector and parameterise it
    s0 = bf.make_beam(matrix.col.random(3, 0.5, 1.5),
                      wavelength=random.uniform(0.8,1.5))
    s0p = BeamParameterisation(s0)

    # apply a random parameter shift
    p_vals = s0p.get_param_vals()
    p_vals = random_param_shift(p_vals, [1000*pi/9, 1000*pi/9, 0.01])
    s0p.set_param_vals(p_vals)

    # compare analytical and finite difference derivatives
    an_ds_dp = s0p.get_ds_dp()
    fd_ds_dp = get_fd_gradients(s0p, [1.e-5 * pi/180, 1.e-5 * pi/180, 1.e-6])

    for j in range(3):
      try:
        assert(approx_equal((fd_ds_dp[j] - an_ds_dp[j]),
                matrix.col((0., 0., 0.)), eps = 1.e-6))
      except Exception:
        failures += 1
        print "for try", i
        print "failure for parameter number", j
        print "with fd_ds_dp = "
        print fd_ds_dp[j]
        print "and an_ds_dp = "
        print an_ds_dp[j]
        print "so that difference fd_ds_dp - an_ds_dp ="
        print fd_ds_dp[j] - an_ds_dp[j]

  if failures == 0: print "OK"
