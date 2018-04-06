from __future__ import absolute_import, division, print_function

import random
import math

import pytest

from dials.algorithms.refinement.parameterisation.goniometer_parameters import \
    GoniometerParameterisation
from dials.algorithms.refinement.refinement_helpers \
    import get_fd_gradients, random_param_shift
from dxtbx.model import Goniometer
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex

def random_gonio():
  # make a random rotation axis with a random setting matrix
  axis = matrix.col(flex.random_double_point_on_sphere())
  fixed_rotation = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
  setting_rotation = matrix.sqr(
      flex.random_double_r3_rotation_matrix())
  goniometer = Goniometer(axis, fixed_rotation, setting_rotation)
  return goniometer

def test():
  goniometer = random_gonio()
  gonp = GoniometerParameterisation(goniometer)

  # Let's do some basic tests. First, can we change parameter values and
  # update the laboratory frame rotation axis?
  e_lab = matrix.col(goniometer.get_rotation_axis())
  gonp.set_param_vals([1000*0.1, 1000*0.1])
  assert matrix.col(goniometer.get_rotation_axis()).angle(e_lab) == pytest.approx(0.1413033)

  # random goniometers and random parameter shifts
  attempts = 1000
  for i in range(attempts):
    # make a random goniometer and parameterise it
    goniometer = random_gonio()
    gonp = GoniometerParameterisation(goniometer)

    # apply a random parameter shift
    p_vals = gonp.get_param_vals()
    p_vals = random_param_shift(p_vals, [1000 * math.pi/9, 1000 * math.pi/9])
    gonp.set_param_vals(p_vals)

    # compare analytical and finite difference derivatives
    an_ds_dp = gonp.get_ds_dp()
    fd_ds_dp = get_fd_gradients(gonp, [1.e-5 * math.pi/180, 1.e-5 * math.pi/180])

    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))
    for j in range(2):
      try:
        assert approx_equal((fd_ds_dp[j] - an_ds_dp[j]),
                null_mat, eps = 1.e-6)
      except Exception:
        print("for try", i)
        print("failure for parameter number", j)
        print("with fd_ds_dp = ")
        print(fd_ds_dp[j])
        print("and an_ds_dp = ")
        print(an_ds_dp[j])
        print("so that difference fd_ds_dp - an_ds_dp =")
        print(fd_ds_dp[j] - an_ds_dp[j])
        raise
