#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function
import logging

logger = logging.getLogger(__name__)

from dials.algorithms.refinement.parameterisation.model_parameters import (
    Parameter,
    ModelParameterisation,
)
from scitbx import matrix
from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge
from dials.algorithms.refinement.refinement_helpers import CrystalOrientationCompose


class CrystalOrientationMixin(object):
    """Mix-in class defining some functionality unique to crystal orientation
    parameterisations that can be shared by static and scan-varying versions"""

    @staticmethod
    def _build_p_list(parameter_type=Parameter):
        """Build the list of parameters, using the parameter_type callback to
        select between versions of the Parameter class"""

        # set up the parameters
        phi1 = parameter_type(0.0, matrix.col((1, 0, 0)), "angle (mrad)", "Phi1")
        phi2 = parameter_type(0.0, matrix.col((0, 1, 0)), "angle (mrad)", "Phi2")
        phi3 = parameter_type(0.0, matrix.col((0, 0, 1)), "angle (mrad)", "Phi3")

        # build the parameter list in a specific,  maintained order
        p_list = [phi1, phi2, phi3]

        return p_list


class CrystalOrientationParameterisation(
    ModelParameterisation, CrystalOrientationMixin
):
    """A parameterisation of the orientation of a Crystal model.

    The Crystal orientation matrix U is parameterised by three Tait-Bryan angles
    expressed in mrad"""

    def __init__(self, crystal, experiment_ids=None):
        """Initialise the CrystalOrientationParameterisation object

        Args:
            crystal: A dxtbx Crystal object to be parameterised.
            experiment_ids (list): The experiment IDs affected by this
                parameterisation. Defaults to None, which is replaced by [0].
        """

        # The state of a crystal orientation parameterisation is an orientation
        # matrix '[U]'. The initial state is a snapshot of the crystal
        # orientation at the time of initialisation '[U0]'. Future states are
        # composed by rotations around axes of the phi-axis frame by Tait-Bryan
        # angles.
        #
        # [U] = [Phi3][Phi2][Phi1][U0]

        # set up the initial state
        if experiment_ids is None:
            experiment_ids = [0]
        istate = matrix.sqr(crystal.get_U())

        # build the parameter list
        p_list = self._build_p_list()

        # set up the base class
        ModelParameterisation.__init__(
            self, crystal, istate, p_list, experiment_ids=experiment_ids
        )

        # call compose to calculate all the derivatives
        self.compose()

        return

    def compose(self):

        # Extract orientation from the initial state
        U0 = self._initial_state

        # extract parameters from the internal list
        phi1, phi2, phi3 = self._param

        # calculate using the helper class
        coc = CrystalOrientationCompose(
            U0, phi1.value, phi1.axis, phi2.value, phi2.axis, phi3.value, phi3.axis
        )

        # compose new state
        self._model.set_U(coc.U())

        # store derivatives
        self._dstate_dp = [coc.dU_dphi1(), coc.dU_dphi2(), coc.dU_dphi3()]

        return

    def get_state(self):

        # only a single crystal is parameterised here, so no multi_state_elt
        # argument is allowed
        return matrix.sqr(self._model.get_U())


class CrystalUnitCellMixin(object):
    """Mix-in class defining some functionality unique to crystal unit cell
    parameterisations that can be shared by static and scan-varying versions"""

    def _build_p_list(self, crystal, parameter_type=Parameter):
        """Build the list of parameters, using the parameter_type callback to
        select between versions of the Parameter class"""

        # Set up symmetrizing object
        S = symmetrize_reduce_enlarge(crystal.get_space_group())
        S.set_orientation(orientation=crystal.get_B())
        X = S.forward_independent_parameters()

        # Set up the independent parameters, with a change of scale
        p_list = [
            parameter_type(e * 1.0e5, name="g_param_%d" % i) for i, e in enumerate(X)
        ]

        return p_list

    def _compose_core(self, raw_vals):

        # obtain metrical matrix parameters on natural scale
        vals = [v * 1.0e-5 for v in raw_vals]

        # set parameter values in the symmetrizing object and obtain new B
        S = symmetrize_reduce_enlarge(self._model.get_space_group())
        S.set_orientation(orientation=self._model.get_B())
        S.forward_independent_parameters()  # Set Bconverter as side-effect
        try:
            newB = matrix.sqr(S.backward_orientation(vals).reciprocal_matrix())
        except RuntimeError as e:

            # write original error to debug log
            logger.debug("Unable to compose the crystal model")
            logger.debug("Original error message: {}".format(str(e)))
            logger.debug("Failing now.")
            raise RuntimeError(
                "Unable to compose the crystal model. Please check that the "
                "experiments match the indexing of the reflections."
            )

        # returns the independent parameters given the set_orientation() B matrix
        # used here for side effects
        S.forward_independent_parameters()

        # get the derivatives of state wrt metrical matrix parameters on the
        # adjusted sale
        dB_dval = [matrix.sqr(e) * 1.0e-5 for e in S.forward_gradients()]

        return newB, dB_dval


class CrystalUnitCellParameterisation(ModelParameterisation, CrystalUnitCellMixin):
    """A parameterisation for the unit cell of a Crystal model.

    The Crystal reciprocal space orthogonalisation matrix B is parameterised
    using up to 6 metrical matrix elements, rescaled by a multiplicative factor.
    """

    def __init__(self, crystal, experiment_ids=None):
        """Initialise the CrystalUnitCellParameterisation object

        Args:
            crystal: A dxtbx Crystal object to be parameterised.
            experiment_ids (list): The experiment IDs affected by this
                parameterisation. Defaults to None, which is replaced by [0].
        """

        # The state of the unit cell parameterisation is the reciprocal space
        # orthogonalisation matrix 'B'. The initial state is irrelevant for
        # this model, as composition of a new B matrix and derivatives can be
        # done with just the values of 6 unit cell parameters, without
        # defining axial directions (which are selected by choice of the PDB
        # convention). For this reason also, the axes of the
        # parameters are irrelevant and are set here to None.

        ### Set up the initial state
        if experiment_ids is None:
            experiment_ids = [0]
        istate = None

        # build the parameter list
        p_list = self._build_p_list(crystal)

        # set up the base class
        ModelParameterisation.__init__(
            self, crystal, istate, p_list, experiment_ids=experiment_ids
        )

        # call compose to calculate all the derivatives
        self.compose()

        return

    def compose(self):

        # calculate new B and derivatives
        newB, self._dstate_dp = self._compose_core([p.value for p in self._param])

        # Now pass new B to the crystal model
        self._model.set_B(newB)

        return

    def get_state(self):

        # only a single crystal is parameterised here, so no multi_state_elt
        # argument is allowed
        return matrix.sqr(self._model.get_B())

    def set_state_uncertainties(self, var_cov, multi_state_elt=None):

        self._model.set_B_covariance(var_cov)

        return
