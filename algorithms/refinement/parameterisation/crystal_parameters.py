#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from model_parameters import Parameter, ModelParameterisation
from dxtbx.model.crystal import crystal_model # implicit import
from scitbx import matrix
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle
from rstbx.symmetry.constraints.parameter_reduction \
    import symmetrize_reduce_enlarge

class CrystalOrientationParameterisation(ModelParameterisation):
  """Parameterisation for crystal orientation, with angles expressed in
  mrad"""

  def __init__(self, crystal, experiment_ids=[0]):

    # The state of a crystal orientation parameterisation is an orientation
    # matrix '[U]'. The initial state is a snapshot of the crystal
    # orientation at the time of initialisation '[U0]'. Future states are
    # composed by rotations around axes of the phi-axis frame by Tait-Bryan
    # angles.
    #
    # [U] = [Phi3][Phi2][Phi1][U0]

    ### Set up the initial state
    istate = crystal.get_U()

    ### Set up the parameters
    phi1 = Parameter(.0, matrix.col((1, 0, 0)), 'angle (mrad)', 'Phi1')
    phi2 = Parameter(.0, matrix.col((0, 1, 0)), 'angle (mrad)', 'Phi2')
    phi3 = Parameter(.0, matrix.col((0, 0, 1)), 'angle (mrad)', 'Phi3')

    # build the parameter list in a specific,  maintained order
    p_list = [phi1, phi2, phi3]

    # set up the base class
    ModelParameterisation.__init__(self, crystal, istate, p_list,
                                   experiment_ids=experiment_ids)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):
    """calculate state and derivatives"""

    # Extract orientation from the initial state
    U0 = self._initial_state

    # extract parameters from the internal list
    phi1, phi2, phi3 = self._param

    # convert to radians
    phi1rad, phi2rad, phi3rad = (phi1.value / 1000., phi2.value / 1000.,
                                 phi3.value / 1000.)

    # compose rotation matrices and their first order derivatives
    Phi1 = (phi1.axis).axis_and_angle_as_r3_rotation_matrix(
                                                        phi1rad, deg=False)
    dPhi1_dphi1 = dR_from_axis_and_angle(phi1.axis, phi1rad, deg=False)

    Phi2 = (phi2.axis).axis_and_angle_as_r3_rotation_matrix(
                                                        phi2rad, deg=False)
    dPhi2_dphi2 = dR_from_axis_and_angle(phi2.axis, phi2rad, deg=False)

    Phi3 = (phi3.axis).axis_and_angle_as_r3_rotation_matrix(
                                                        phi3rad, deg=False)
    dPhi3_dphi3 = dR_from_axis_and_angle(phi3.axis, phi3rad, deg=False)

    Phi21 = Phi2 * Phi1
    Phi321 = Phi3 * Phi21

    ### Compose new state

    newU = Phi321 * U0
    self._model.set_U(newU)

    ### calculate derivatives of the state wrt parameters
    dU_dphi1 = Phi3 * Phi2 * dPhi1_dphi1 * U0
    dU_dphi2 = Phi3 * dPhi2_dphi2 * Phi1 * U0
    dU_dphi3 = dPhi3_dphi3 * Phi21 * U0

    # convert to mrad and store
    self._dstate_dp = [dU_dphi1 / 1000., dU_dphi2 / 1000., dU_dphi3 / 1000.]

    return

  def get_state(self):

    # only a single crystal is parameterised here, so no multi_state_elt
    # argument is allowed
    return matrix.sqr(self._model.get_U())

class CrystalUnitCellParameterisation(ModelParameterisation):
  """Parameterisation for the unit cell"""

  def __init__(self, crystal, experiment_ids=[0]):

    # The state of the unit cell parameterisation is the reciprocal space
    # orthogonalisation matrix 'B'. The initial state is irrelevant for
    # this model, as composition of a new B matrix and derivatives can be
    # done with just the values of 6 unit cell parameters, without
    # defining axial directions (which are selected by choice of the PDB
    # convention). For this reason also, the axes of the
    # parameters are irrelevant and are set here to None.

    ### Set up the initial state
    istate = None

    ### Set up symmetrizing object
    self._S = symmetrize_reduce_enlarge(crystal.get_space_group())
    self._S.set_orientation(orientation=crystal.get_B())
    X = self._S.forward_independent_parameters()
    dB_dp = self._S.forward_gradients()
    B = self._S.backward_orientation(independent=X).reciprocal_matrix()

    ### Set up the independent parameters, with a change of scale
    p_list = [Parameter(e * 1.e5, name = "g_param_%d" % i) \
              for i, e in enumerate(X)]

    # set up the base class
    ModelParameterisation.__init__(self, crystal, istate, p_list,
                                   experiment_ids=experiment_ids)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):
    """calculate state and derivatives"""

    # obtain parameters on natural scale
    p_vals = [p.value / 1.e5 for p in self._param]

    # set parameter values in the symmetrizing object and obtain new B
    newB = matrix.sqr(
            self._S.backward_orientation(p_vals).reciprocal_matrix())

    # Now pass new B to the crystal model
    self._model.set_B(newB)

    # returns the independent parameters given the set_orientation() B
    # matrix. Used here for side effects
    self._S.forward_independent_parameters()

    # get the gradients on the adjusted scale
    self._dstate_dp = [matrix.sqr(e) / 1.e5 \
                       for e in self._S.forward_gradients()]

    return

  def get_state(self):

    # only a single crystal is parameterised here, so no multi_state_elt
    # argument is allowed
    return matrix.sqr(self._model.get_B())

  def set_state_uncertainties(self, var_cov, multi_state_elt=None):

    from dials.algorithms.refinement.refinement_helpers import \
        matrix_inverse_error_propagation

    # var_cov is the covariance matrix of elements of the B matrix. We
    # need to construct the covariance matrix of elements of the
    # transpose of B. The vector of elements of B is related to the
    # vector of elements of its transpose by a permutation, P.
    P = matrix.sqr((1,0,0,0,0,0,0,0,0,
                    0,0,0,1,0,0,0,0,0,
                    0,0,0,0,0,0,1,0,0,
                    0,1,0,0,0,0,0,0,0,
                    0,0,0,0,1,0,0,0,0,
                    0,0,0,0,0,0,0,1,0,
                    0,0,1,0,0,0,0,0,0,
                    0,0,0,0,0,1,0,0,0,
                    0,0,0,0,0,0,0,0,1))

    # We can use P to replace var_cov with the covariance matrix of
    # elements of the transpose of B.
    var_cov = P*var_cov*P

    # From B = (O^-1)^T we can convert this
    # to the covariance matrix of the real space orthogonalisation matrix
    B = self.get_state()
    Bt = B.transpose()
    O = Bt.inverse()
    cov_O = matrix_inverse_error_propagation(Bt, var_cov)

    # The real space unit cell vectors are given by
    vec_a = (O * matrix.col((1,0,0)))
    vec_b = (O * matrix.col((0,1,0)))
    vec_c = (O * matrix.col((0,0,1)))

    a, b, c = vec_a.length(), vec_b.length(), vec_c.length()

    # The estimated errors are calculated by error propagation from cov_O
    jacobian = matrix.rec(
      (vec_a[0]/a, 0, 0, vec_a[1]/a, 0, 0, vec_a[2]/a, 0, 0,
       0, vec_b[0]/b, 0, 0, vec_b[1]/b, 0, 0, vec_b[2]/b, 0,
       0, 0, vec_c[0]/c, 0, 0, vec_c[1]/c, 0, 0, vec_c[2]/c), (3, 9))
    jacobian_t = jacobian.transpose()

    # this is a 3*3 covariance matrix, we want the variances, down the diagonal
    cov_cell = (jacobian * cov_O * jacobian_t)
    var_a, var_b, var_c = cov_cell[0], cov_cell[4], cov_cell[8]

    # FIXME These estimates are not tested (how can we even do so? Use semi-
    # synthetic datasets perhaps and look at true scatter?) plus we have nowhere
    # to store this information. So for now just output to the debug log
    from logging import debug
    from math import sqrt
    debug("Refined cell lengths and estimated standard deviations:")
    debug("a: {0:f} +/- ({1:f})".format(a, sqrt(var_a)))
    debug("b: {0:f} +/- ({1:f})".format(b, sqrt(var_b)))
    debug("c: {0:f} +/- ({1:f})".format(c, sqrt(var_c)))

    # TODO cell angles. How?

    return
