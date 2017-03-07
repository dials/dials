#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
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
from dials.algorithms.refinement.parameterisation.crystal_parameters \
    import CrystalOrientationMixin, CrystalUnitCellMixin
from dials.algorithms.refinement.refinement_helpers import CrystalOrientationCompose

class ScanVaryingCrystalOrientationParameterisation(
  ScanVaryingModelParameterisation, CrystalOrientationMixin):
  """Scan-varying parameterisation for crystal orientation, with angles
  expressed in mrad"""

  def __init__(self, crystal, t_range, num_intervals, experiment_ids=None):
    if experiment_ids is None:
      experiment_ids = [0]

    # The state of a scan varying crystal orientation parameterisation
    # is an orientation
    # matrix '[U](t)', expressed as a function of image number 't'
    # in a sequential scan.
    #
    # The initial state is a snapshot of the crystal orientation
    # at the point of initialisation '[U0]', which is independent of
    # image number.
    #
    # Future states are composed by
    # rotations around axes of the phi-axis frame by Tait-Bryan angles.
    #
    # [U](t) = [Phi3](t)[Phi2](t)[Phi1](t)[U0]

    # Set up the smoother
    smoother = GaussianSmoother(t_range, num_intervals)
    nv = smoother.num_values()

    # Set up the initial state
    istate = matrix.sqr(crystal.get_U())
    self._U_at_t = istate

    # Factory function to provide to _build_p_list
    def parameter_type(value, axis, ptype, name):
      return ScanVaryingParameterSet(value, nv, axis, ptype, name)

    # Build the parameter list
    p_list = self._build_p_list(parameter_type)

    # Set up the base class
    ScanVaryingModelParameterisation.__init__(self, crystal, istate,
                                              p_list, smoother,
                                              experiment_ids=experiment_ids)

    return

  def compose(self, t):
    """calculate state and derivatives for model at image number t"""

    # Extract orientation from the initial state
    U0 = self._initial_state

    # extract parameter sets from the internal list
    phi1_set, phi2_set, phi3_set = self._param

    # extract angles and other data at time t using the smoother
    phi1, phi1_weights, phi1_sumweights = self._smoother.value_weight(t, phi1_set)
    phi2, phi2_weights, phi2_sumweights = self._smoother.value_weight(t, phi2_set)
    phi3, phi3_weights, phi3_sumweights = self._smoother.value_weight(t, phi3_set)

    # calculate derivatives of angles wrt underlying parameters.
    # FIXME write up notes in orange notebook
    dphi1_dp = phi1_weights * (1. / phi1_sumweights)
    dphi2_dp = phi2_weights * (1. / phi2_sumweights)
    dphi3_dp = phi3_weights * (1. / phi3_sumweights)

    # calculate state and derivatives using the helper class
    coc = CrystalOrientationCompose(U0, phi1, phi1_set.axis,
                                    phi2, phi2_set.axis,
                                    phi3, phi3_set.axis)
    self._U_at_t = coc.U()
    dU_dphi1 = coc.dU_dphi1()
    dU_dphi2 = coc.dU_dphi2()
    dU_dphi3 = coc.dU_dphi3()

    # calculate derivatives of state wrt underlying parameters
    dU_dp1 = [None] * dphi1_dp.size
    for (i, v) in dphi1_dp: dU_dp1[i] = dU_dphi1 * v
    dU_dp2 = [None] * dphi2_dp.size
    for (i, v) in dphi2_dp: dU_dp2[i] = dU_dphi2 * v
    dU_dp3 = [None] * dphi3_dp.size
    for (i, v) in dphi3_dp: dU_dp3[i] = dU_dphi3 * v

    # store derivatives as list-of-lists
    self._dstate_dp = [dU_dp1, dU_dp2, dU_dp3]

    return

  def get_state(self):
    """Return crystal orientation matrix [U] at image number t"""

    # only a single crystal is parameterised here, so no multi_state_elt
    # argument is allowed

    return self._U_at_t


class ScanVaryingCrystalUnitCellParameterisation(
  ScanVaryingModelParameterisation, CrystalUnitCellMixin):
  """Scan-varying parameterisation for the crystal unit cell"""

  def __init__(self, crystal, t_range, num_intervals, experiment_ids=None):
    from scitbx import matrix
    if experiment_ids is None:
      experiment_ids = [0]

    # The state of a scan-varying unit cell parameterisation is the
    # reciprocal space orthogonalisation matrix '[B](t)', expressed as a
    # function of image number 't' in a sequential scan.

    # Other comments from CrystalUnitCellParameterisation are relevant here

    # Set up the smoother
    smoother = GaussianSmoother(t_range, num_intervals)
    nv = smoother.num_values()

    # Set up the initial state
    istate = None
    self._B_at_t = matrix.sqr(crystal.get_B())

    # Factory function to provide to _build_p_list
    def parameter_type(value, name):
      return ScanVaryingParameterSet(value, nv, name=name)

    # Build the parameter list
    p_list = self._build_p_list(crystal, parameter_type)

    # Set up the base class
    ScanVaryingModelParameterisation.__init__(self, crystal, istate,
                                              p_list, smoother,
                                              experiment_ids=experiment_ids)

    return

  def compose(self, t):
    """calculate state and derivatives for model at image number t"""

    # extract values and weights at time t using the smoother
    vals, weights, sumweights = zip(*(self._smoother.value_weight(t,
      pset) for pset in self._param))

    # calculate derivatives of metrical matrix parameters wrt underlying
    # scan-varying parameters
    inv_sumw = [1. / sw for sw in sumweights]
    dvals_dp =  [e * isw for e, isw in zip(weights, inv_sumw)]

    # calculate new B and derivatives
    self._B_at_t, dB_dval = self._compose_core(vals)

    # calculate derivatives of state wrt underlying parameters
    self._dstate_dp = [[b * e for e in a.as_dense_vector()] for a, b in zip(dvals_dp, dB_dval)]
    self._dstate_dp = [[None] * e.size for e in dvals_dp]
    for i, (dv, dB) in enumerate(zip(dvals_dp, dB_dval)):
      for j, e in dv: self._dstate_dp[i][j] = e * dB

    return

  def get_state(self):
    """Return crystal orthogonalisation matrix [B] at image number t"""

    # only a single crystal is parameterised here, so no multi_state_elt
    # argument is allowed

    return self._B_at_t
