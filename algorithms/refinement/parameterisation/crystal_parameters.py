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
from dials.model.experiment.crystal_model import Crystal
from scitbx import matrix
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle
from rstbx.symmetry.constraints.parameter_reduction \
    import symmetrize_reduce_enlarge

class CrystalOrientationParameterisation(ModelParameterisation):
  """Parameterisation for crystal orientation, with angles expressed in
  mrad"""

  def __init__(self, crystal):

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

    # set up the list of model objects being parameterised (here
    # just a single crystal model)
    models = [crystal]

    # set up the base class
    ModelParameterisation.__init__(self, models, istate, p_list)

    # call compose to calculate all the derivatives
    self.compose()

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
    self._models[0].set_U(newU)

    ### calculate derivatives of the state wrt parameters
    dU_dphi1 = Phi3 * Phi2 * dPhi1_dphi1 * U0
    dU_dphi2 = Phi3 * dPhi2_dphi2 * Phi1 * U0
    dU_dphi3 = dPhi3_dphi3 * Phi21 * U0

    # convert to mrad and store
    self._dstate_dp = [dU_dphi1 / 1000., dU_dphi2 / 1000., dU_dphi3 / 1000.]

    return

  def get_state(self):
    return matrix.sqr(self._models[0].get_U())

class CrystalUnitCellParameterisation(ModelParameterisation):
  """Parameterisation for the unit cell"""

  def __init__(self, crystal):

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

    # set up the list of model objects being parameterised (here
    # just a single crystal model)
    models = [crystal]

    # set up the base class
    ModelParameterisation.__init__(self, models, istate, p_list)

    # call compose to calculate all the derivatives
    self.compose()

  def compose(self):
    """calculate state and derivatives"""

    # obtain parameters on natural scale
    p_vals = [p.value / 1.e5 for p in self._param]

    # set parameter values in the symmetrizing object and obtain new B
    newB = matrix.sqr(
            self._S.backward_orientation(p_vals).reciprocal_matrix())

    # Now pass new B to the crystal model
    self._models[0].set_B(newB)

    # returns the independent parameters given the set_orientation() B
    # matrix. Used here for side effects
    self._S.forward_independent_parameters()

    # get the gradients on the adjusted scale
    self._dstate_dp = [matrix.sqr(e) / 1.e5 \
                       for e in self._S.forward_gradients()]

    return

  def get_state(self):
    return matrix.sqr(self._models[0].get_B())
