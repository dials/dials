#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
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
from dials.algorithms.refinement.parameterisation.beam_parameters import BeamMixin

class ScanVaryingBeamParameterisation(
  ScanVaryingModelParameterisation, BeamMixin):
  """A scan-varying parameterisation for the beam with angles expressed in
  mrad and wavenumber in inverse Angstroms."""

  def __init__(self, beam, t_range, num_intervals, goniometer=None,
      experiment_ids=None):

    if experiment_ids is None:
      experiment_ids = [0]

    # The state of a scan varying beam parameterisation is a vector
    # '[s0](t)', expressed as a function of image number 't'
    # in a sequential scan.
    #
    # The initial state is a snapshot of the beam unit direction
    # at the point of initialisation '[s0dir]', which is independent of
    # image number.
    #
    # Future states are composed by rotations around axes orthogonal to the
    # initial direction and scaling by the wavenumber.
    #
    # [s0](t) = nu(t) * [Mu2](t)[Mu1](t)[s0dir]

    # Set up the smoother
    smoother = GaussianSmoother(t_range, num_intervals)
    nv = smoother.num_values()

    # Set up the initial state
    s0 = matrix.col(beam.get_s0())
    s0dir = matrix.col(beam.get_unit_s0())
    istate = s0dir
    self._s0_at_t = s0

    # Factory function to provide to _build_p_list
    def parameter_type(value, axis, ptype, name):
      return ScanVaryingParameterSet(value, nv, axis, ptype, name)

    # Build the parameter list
    p_list = self._build_p_list(s0, goniometer, parameter_type=parameter_type)

    # Set up the base class
    ScanVaryingModelParameterisation.__init__(self, beam, istate,
                                              p_list, smoother,
                                              experiment_ids=experiment_ids)

    return

  def compose(self, t):
    """calculate state and derivatives for model at image number t"""

    # extract direction from the initial state
    is0 = self._initial_state

    # extract parameter sets from the internal list
    mu1_set, mu2_set, nu_set = self._param

    # extract data at time t using the smoother
    mu1, mu1_weights, mu1_sumweights = self._smoother.value_weight(t, mu1_set)
    mu2, mu2_weights, mu2_sumweights = self._smoother.value_weight(t, mu2_set)
    nu, nu_weights, nu_sumweights = self._smoother.value_weight(t, nu_set)

    # calculate derivatives of values at image t wrt underlying parameters.
    dmu1_dp = mu1_weights * (1. / mu1_sumweights)
    dmu2_dp = mu2_weights * (1. / mu2_sumweights)
    dnu_dp = nu_weights * (1. / nu_sumweights)

    # calculate new s0 and its derivatives wrt the values mu1, mu2 and nu
    self._s0_at_t, ds0_dval = self._compose_core(is0, mu1, mu2, nu,
      mu1_axis=mu1_set.axis, mu2_axis=mu2_set.axis)

    # calculate derivatives of state wrt underlying smoother parameters
    ds0_dp1 = [None] * dmu1_dp.size
    for (i, v) in dmu1_dp: ds0_dp1[i] = ds0_dval[0] * v
    ds0_dp2 = [None] * dmu2_dp.size
    for (i, v) in dmu2_dp: ds0_dp2[i] = ds0_dval[1] * v
    ds0_dp3 = [None] * dnu_dp.size
    for (i, v) in dnu_dp: ds0_dp3[i] = ds0_dval[2] * v

    # store derivatives as list-of-lists
    self._dstate_dp = [ds0_dp1, ds0_dp2, ds0_dp3]

    return

  def get_state(self):
    """Return beam vector [s0] at image number t"""

    # only a single beam is parameterised here, so no multi_state_elt
    # argument is allowed

    return self._s0_at_t
