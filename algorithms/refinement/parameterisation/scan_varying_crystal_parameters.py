#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from logging import debug
from scitbx import matrix
from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
        import ScanVaryingParameterSet, \
               ScanVaryingModelParameterisation, \
               GaussianSmoother
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle

class ScanVaryingCrystalOrientationParameterisation(ScanVaryingModelParameterisation):
  """A work-in-progress time-dependent parameterisation for crystal
  orientation, with angles expressed in mrad"""

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
    istate = crystal.get_U()

    # Set up the parameters
    phi1 = ScanVaryingParameterSet(0.0, nv,
                        matrix.col((1., 0., 0.)), 'angle (mrad)', 'Phi1')
    phi2 = ScanVaryingParameterSet(0.0, nv,
                        matrix.col((0., 1., 0.)), 'angle (mrad)', 'Phi2')
    phi3 = ScanVaryingParameterSet(0.0, nv,
                        matrix.col((0., 0., 1.)), 'angle (mrad)', 'Phi3')

    # Build the list of parameter sets in a specific, maintained order
    p_list = [phi1, phi2, phi3]

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
    dphi1_dp = [e / phi1_sumweights for e in phi1_weights]
    dphi2_dp = [e / phi2_sumweights for e in phi2_weights]
    dphi3_dp = [e / phi3_sumweights for e in phi3_weights]

    # convert angles to radians
    phi1rad, phi2rad, phi3rad = (phi1 / 1000., phi2 / 1000.,
                                 phi3 / 1000.)

    # compose rotation matrices and their first order derivatives wrt angle
    Phi1 = (phi1_set.axis).axis_and_angle_as_r3_rotation_matrix(phi1rad, deg=False)
    dPhi1_dphi1 = dR_from_axis_and_angle(phi1_set.axis, phi1rad, deg=False)

    Phi2 = (phi2_set.axis).axis_and_angle_as_r3_rotation_matrix(phi2rad, deg=False)
    dPhi2_dphi2 = dR_from_axis_and_angle(phi2_set.axis, phi2rad, deg=False)

    Phi3 = (phi3_set.axis).axis_and_angle_as_r3_rotation_matrix(phi3rad, deg=False)
    dPhi3_dphi3 = dR_from_axis_and_angle(phi3_set.axis, phi3rad, deg=False)

    Phi21 = Phi2 * Phi1
    Phi321 = Phi3 * Phi21

    # Compose new state
    self._U_at_t = Phi321 * U0

    # calculate derivatives of the state wrt angle, convert back to mrad
    dU_dphi1 = Phi3 * Phi2 * dPhi1_dphi1 * U0 / 1000.
    dU_dphi2 = Phi3 * dPhi2_dphi2 * Phi1 * U0 / 1000.
    dU_dphi3 = dPhi3_dphi3 * Phi21 * U0 / 1000.

    # calculate derivatives of state wrt underlying parameters
    dU_dp1 = [dU_dphi1 * e for e in dphi1_dp]
    dU_dp2 = [dU_dphi2 * e for e in dphi2_dp]
    dU_dp3 = [dU_dphi3 * e for e in dphi3_dp]

    # store derivatives as list-of-lists
    self._dstate_dp = [dU_dp1, dU_dp2, dU_dp3]

    return

  def get_state(self):
    """Return crystal orientation matrix [U] at image number t"""

    # only a single crystal is parameterised here, so no multi_state_elt
    # argument is allowed

    return self._U_at_t


class ScanVaryingCrystalUnitCellParameterisation(ScanVaryingModelParameterisation):
  """A work-in-progress time-dependent parameterisation for the crystal
  unit cell"""

  def __init__(self, crystal, t_range, num_intervals, experiment_ids=None):
    if experiment_ids is None:
      experiment_ids = [0]

    # The state of a scan-varying unit cell parameterisation is the
    # reciprocal space orthogonalisation matrix '[B](t)', expressed as a
    # function of image number 't' in a sequential scan.

    # Other comments from CrystalUnitCellParameterisation are relevant here

    # Set up the smoother
    smoother = GaussianSmoother(t_range, num_intervals)
    nv = smoother.num_values()

    ### Set up the initial state
    istate = None

    ### Set up symmetrizing object
    self._S = symmetrize_reduce_enlarge(crystal.get_space_group())
    self._S.set_orientation(orientation=crystal.get_B())
    X = self._S.forward_independent_parameters()
    dB_dp = self._S.forward_gradients()
    B = self._S.backward_orientation(independent=X).reciprocal_matrix()

    ### Set up the independent parameters, with a change of scale
    p_list = [ScanVaryingParameterSet(e * 1.e5, nv, name = "g_param_%d" % i) \
              for i, e in enumerate(X)]

    # Set up the base class
    ScanVaryingModelParameterisation.__init__(self, crystal, istate,
                                              p_list, smoother,
                                              experiment_ids=experiment_ids)

    return

  def compose(self, t):
    """calculate state and derivatives for model at image number t"""

    # extract values and weights at time t using the smoother
    data = [self._smoother.value_weight(t, pset) for pset in self._param]

    # obtain metrical matrix parameters on natural scale
    vals = [val / 1.e5 for val, weights, sumweight in data]

    # calculate derivatives of metrical matrix parameters wrt underlying
    # scan-varying parameters
    dvals_dp =  [tuple([e / sw for e in w]) for v, w, sw in data]

    # set parameter values in the symmetrizing object and obtain new B
    try:
      self._B_at_t = matrix.sqr(
                  self._S.backward_orientation(vals).reciprocal_matrix())
    except RuntimeError as e:
      from libtbx.utils import Sorry
      # write original error to debug log
      debug('Unable to compose the crystal model')
      debug('Original error message: {0}'.format(str(e)))
      debug('Failing now.')
      raise Sorry('Unable to compose the crystal model. Please check that the '
                  'experiments match the indexing of the reflections.')

    # returns the independent parameters given the set_orientation() B matrix
    # used here for side effects
    self._S.forward_independent_parameters()

    # get the derivatives of state wrt metrical matrix parameters on the
    # adjusted scale
    dB_dval = [matrix.sqr(e) / 1.e5 \
                       for e in self._S.forward_gradients()]

    # calculate derivatives of state wrt underlying parameters
    self._dstate_dp = [[b * e for e in a] for a, b in zip(dvals_dp, dB_dval)]

    return

  def get_state(self):
    """Return crystal orthogonalisation matrix [B] at image number t"""

    # only a single crystal is parameterised here, so no multi_state_elt
    # argument is allowed

    return self._B_at_t
