#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from collections import namedtuple
from scitbx import matrix
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
        import ScanVaryingParameterSet, \
               ScanVaryingModelParameterisation, \
               GaussianSmoother
from dials.algorithms.refinement.parameterisation.detector_parameters import \
  DetectorMixin

# bucket to provide value and axis attributes to model the Parameter class
Param = namedtuple('Param', ['value', 'axis'])

class ScanVaryingDetectorParameterisationSinglePanel(
  ScanVaryingModelParameterisation, DetectorMixin):
  """A scan-varying parameterisation for a single abstract panel plane, with
  angles expressed in mrad"""

  def __init__(self, detector, t_range, num_intervals, experiment_ids=None):

    if experiment_ids is None:
      experiment_ids = [0]

    # Set up the smoother
    smoother = GaussianSmoother(t_range, num_intervals)
    nv = smoother.num_values()

    # Factory function to provide to _init_core
    def parameter_type(value, axis, ptype, name):
      return ScanVaryingParameterSet(value, nv, axis, ptype, name)

    # Set up the initial state and parameter list
    dat = self._init_core(detector, parameter_type)

    self._d_at_t = matrix.sqr(detector[0].get_d_matrix())

    # set up the base class
    ScanVaryingModelParameterisation.__init__(self, detector, dat['istate'],
      dat['p_list'], smoother, experiment_ids=experiment_ids)

    return

  def compose(self, t):
    """calculate state and derivatives for model at image number t"""

    # extract parameter sets from the internal list
    dist_set, shift1_set, shift2_set, tau1_set, tau2_set, tau3_set = self._param

    # extract data at time t using the smoother
    dist, dist_weights, dist_sumweights = self._smoother.value_weight(
      t, dist_set)
    shift1, shift1_weights, shift1_sumweights = self._smoother.value_weight(
      t, shift1_set)
    shift2, shift2_weights, shift2_sumweights = self._smoother.value_weight(
      t, shift2_set)
    tau1, tau1_weights, tau1_sumweights = self._smoother.value_weight(
      t, tau1_set)
    tau2, tau2_weights, tau2_sumweights = self._smoother.value_weight(
      t, tau2_set)
    tau3, tau3_weights, tau3_sumweights = self._smoother.value_weight(
      t, tau3_set)

    # calculate derivatives of values at image t wrt underlying parameters.
    ddist_dp = dist_weights * (1. / dist_sumweights)
    dshift1_dp = shift1_weights * (1. / shift1_sumweights)
    dshift2_dp = shift2_weights * (1. / shift2_sumweights)
    dtau1_dp = tau1_weights * (1. / tau1_sumweights)
    dtau2_dp = tau2_weights * (1. / tau2_sumweights)
    dtau3_dp = tau3_weights * (1. / tau3_sumweights)

    # calculate new [d] and derivatives wrt the values dist, shift1, shift2,
    # tau1, tau2 and tau3
    new_state, dd_dval = self._compose_core(
      Param(dist, dist_set.axis),
      Param(shift1, shift1_set.axis),
      Param(shift2, shift2_set.axis),
      Param(tau1, tau1_set.axis),
      Param(tau2, tau2_set.axis),
      Param(tau3, tau3_set.axis))

    # 'fast' axis elements
    f1, f2, f3 = new_state['d1']

    # 'slow' axis elements
    s1, s2, s3 = new_state['d2']

    # origin vector elements
    o1, o2, o3 = new_state['origin']

    self._d_at_t = matrix.sqr((f1, s1, o1,
                               f2, s2, o2,
                               f3, s3, o3))

    # calculate derivatives of state wrt underlying smoother parameters
    dd_dp1 = [None] * ddist_dp.size
    for (i, v) in ddist_dp: dd_dp1[i] = dd_dval[0] * v
    dd_dp2 = [None] * dshift1_dp.size
    for (i, v) in dshift1_dp: dd_dp2[i] = dd_dval[1] * v
    dd_dp3 = [None] * dshift2_dp.size
    for (i, v) in dshift2_dp: dd_dp3[i] = dd_dval[2] * v
    dd_dp4 = [None] * dtau1_dp.size
    for (i, v) in dtau1_dp: dd_dp4[i] = dd_dval[3] * v
    dd_dp5 = [None] * dtau2_dp.size
    for (i, v) in dtau2_dp: dd_dp5[i] = dd_dval[4] * v
    dd_dp6 = [None] * dtau3_dp.size
    for (i, v) in dtau3_dp: dd_dp6[i] = dd_dval[5] * v

    # store derivatives as list-of-lists
    self._dstate_dp = [dd_dp1, dd_dp2, dd_dp3, dd_dp4, dd_dp5, dd_dp6]

    return

  def get_state(self):
    """Return detector matrix [d] at image number t"""
    # only a single panel exists, so no multi_state_elt argument is allowed
    return self._d_at_t
