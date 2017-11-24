#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division
from scitbx import matrix
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
    import ScanVaryingParameterSet, \
           ScanVaryingModelParameterisation, \
           GaussianSmoother
from dials.algorithms.refinement.parameterisation.goniometer_parameters \
    import GoniometerMixin

class ScanVaryingGoniometerParameterisation(
  ScanVaryingModelParameterisation, GoniometerMixin):
  """A scan-varying parameterisation for the setting rotation of a goniometer
  with angles expressed in mrad."""

  def __init__(self, goniometer, t_range, num_intervals, beam=None,
      experiment_ids=None):

    if experiment_ids is None:
      experiment_ids = [0]

    # The state of a scan varying goniometer parameterisation is a matrix
    # '[S](t)', expressed as a function of image number 't'
    # in a sequential scan.
    #
    # The initial state is a snapshot of the setting matrix
    # at the point of initialisation '[S0]', which is independent of
    # image number.
    #
    # Future states are composed by rotations around two axes orthogonal to the
    # initial spindle axis direction.
    #
    # [S](t) = [G2](t)[G1](t)[S0]

    # Set up the smoother
    smoother = GaussianSmoother(t_range, num_intervals)
    nv = smoother.num_values()

    # Set up the initial state
    e_lab = matrix.col(goniometer.get_rotation_axis())
    istate = matrix.sqr(goniometer.get_setting_rotation())
    self._S_at_t = istate

    # Factory function to provide to _build_p_list
    def parameter_type(value, axis, ptype, name):
      return ScanVaryingParameterSet(value, nv, axis, ptype, name)

    # Build the parameter list
    p_list = self._build_p_list(e_lab, beam, parameter_type=parameter_type)

    # Set up the base class
    ScanVaryingModelParameterisation.__init__(self, goniometer, istate,
                                              p_list, smoother,
                                              experiment_ids=experiment_ids)

    return

  def compose(self, t):
    """calculate state and derivatives for model at image number t"""

    # Extract setting matrix from the initial state
    iS0 = self._initial_state

    # extract parameter sets from the internal list
    gamma1_set, gamma2_set = self._param

    # extract angles and other data at time t using the smoother
    gamma1, gamma1_weights, gamma1_sumweights = self._smoother.value_weight(t, gamma1_set)
    gamma2, gamma2_weights, gamma2_sumweights = self._smoother.value_weight(t, gamma2_set)

    # calculate derivatives of angles wrt underlying parameters.
    # FIXME write up notes in orange notebook
    dgamma1_dp = gamma1_weights * (1. / gamma1_sumweights)
    dgamma2_dp = gamma2_weights * (1. / gamma2_sumweights)

    self._S_at_t, dS_dval = self._compose_core(iS0, gamma1, gamma2,
      gamma1_axis=gamma1_set.axis, gamma2_axis=gamma2_set.axis)
    print self._S_at_t, t

    # calculate derivatives of state wrt underlying smoother parameters
    dS_dp1 = [None] * dgamma1_dp.size
    for (i, v) in dgamma1_dp: dS_dp1[i] = dS_dval[0] * v
    dS_dp2 = [None] * dgamma2_dp.size
    for (i, v) in dgamma2_dp: dS_dp2[i] = dS_dval[1] * v

    # store derivatives as list-of-lists
    self._dstate_dp = [dS_dp1, dS_dp2]

    return

  def get_state(self):
    """Return setting matrix [S] at image number t"""

    # only a single goniometer is parameterised here, so no multi_state_elt
    # argument is allowed

    return self._S_at_t
