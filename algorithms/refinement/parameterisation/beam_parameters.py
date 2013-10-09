#!/usr/bin/env python

from __future__ import division
from model_parameters import Parameter, ModelParameterisation
from scitbx import matrix
from dials.algorithms.refinement import dR_from_axis_and_angle

class BeamParameterisationOrientation(ModelParameterisation):
    '''Implementation of parameterisation for the beam (direction only)
    with angles expressed in mrad.

    Pass in a goniometer (if present) to ensure consistent definition of the
    beam rotation angles with respect to the spindle-beam plane.'''

    def __init__(self, beam, goniometer = None):

        # The state of the beam model consists of the orientation of the
        # s0 vector that it is modelling. The initial state is a snapshot
        # of the s0 vector at the point of initialisation. Future states are
        # composed by rotations around axes perpendicular to that direction.
        #
        # The 'models' attribute refers to the beam vector contained by this
        # model.

        ### Set up the initial state
        s0 = matrix.col(beam.get_s0())
        s0dir = matrix.col(beam.get_unit_s0())
        istate = s0dir

        ### Set up the parameters
        if goniometer:
            spindle = matrix.col(goniometer.get_rotation_axis())
            s0_plane_dir2 = s0.cross(spindle).normalize()
            s0_plane_dir1 = s0_plane_dir2.cross(s0).normalize()
        else:
            s0_plane_dir1 = s0.ortho().normalize()
            s0_plane_dir2 = s0.cross(s0_plane_dir1).normalize()

        # rotation around s0_plane_dir1
        mu1 = Parameter(.0, s0_plane_dir1, 'angle (mrad)', 'Mu1')
        # rotation around s0_plane_dir2
        mu2 = Parameter(.0, s0_plane_dir2, 'angle (mrad)', 'Mu2')

        # build the parameter list in a specific,  maintained order
        p_list = [mu1, mu2]

        # set up the list of model objects being parameterised (here
        # just a single beam model)
        models = [beam]

        # set up the base class
        ModelParameterisation.__init__(self, models, istate, p_list)

        # call compose to calculate all the derivatives
        self.compose()

    def compose(self):

        # extract direction from the initial state
        is0 = self._initial_state

        # extract parameters from the internal list
        mu1, mu2 = self._param

        # convert to radians
        mu1rad, mu2rad = mu1.value / 1000., mu2.value / 1000.

        # compose rotation matrices and their first order derivatives
        Mu1 = (mu1.axis).axis_and_angle_as_r3_rotation_matrix(mu1rad, deg=False)
        dMu1_dmu1 = dR_from_axis_and_angle(mu1.axis, mu1rad, deg=False)

        Mu2 = (mu2.axis).axis_and_angle_as_r3_rotation_matrix(mu2rad, deg=False)
        dMu2_dmu2 = dR_from_axis_and_angle(mu2.axis, mu2rad, deg=False)

        Mu21 = Mu2 * Mu1

        ### Compose new state
        s0_new_dir = (Mu21 * is0).normalize()

        # now update the model with its new orientation
        s0mag = 1. / self._models[0].get_wavelength()
        self._models[0].set_unit_s0(s0mag * s0_new_dir)

        ### calculate derivatives of the state wrt parameters
        # derivative wrt mu1
        dMu21_dmu1 = Mu2 * dMu1_dmu1
        ds0_new_dir_dmu1 = dMu21_dmu1 * is0

        # derivative wrt mu2
        dMu21_dmu2 = dMu2_dmu2 * Mu1
        ds0_new_dir_dmu2 = dMu21_dmu2 * is0

        ### calculate derivatives of the attached beam vector, converting
        ### parameters back to mrad, and store
        # derivative wrt mu1
        self._dstate_dp[0] = ds0_new_dir_dmu1 * s0mag / 1000.
        # derivative wrt mu2
        self._dstate_dp[1] = ds0_new_dir_dmu2 * s0mag / 1000.

        return

    def get_state(self):
        return matrix.col(self._models[0].get_s0())
