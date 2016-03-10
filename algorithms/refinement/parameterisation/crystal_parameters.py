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
from logging import debug
from math import sqrt, pi, acos
from model_parameters import Parameter, ModelParameterisation
from dxtbx.model.crystal import crystal_model # implicit import
from scitbx import matrix
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle
from rstbx.symmetry.constraints.parameter_reduction \
    import symmetrize_reduce_enlarge
from dials.algorithms.refinement.refinement_helpers \
  import CrystalOrientationCompose

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

    # calculate using the helper class
    coc = CrystalOrientationCompose(U0, phi1.value, phi1.axis,
                                    phi2.value, phi2.axis,
                                    phi3.value, phi3.axis)

    # compose new state
    self._model.set_U(coc.U())

    # store derivatives
    self._dstate_dp = [coc.dU_dphi1(), coc.dU_dphi2(), coc.dU_dphi3()]

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
    try:
      newB = matrix.sqr(
            self._S.backward_orientation(p_vals).reciprocal_matrix())
    except RuntimeError as e:
      from libtbx.utils import Sorry
      # write original error to debug log
      debug('Unable to compose the crystal model')
      debug('Original error message: {0}'.format(str(e)))
      debug('Failing now.')
      raise Sorry('Unable to compose the crystal model. Please check that the '
                  'experiments match the indexing of the reflections.')

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

    self._model.set_B_covariance(var_cov)

    return
