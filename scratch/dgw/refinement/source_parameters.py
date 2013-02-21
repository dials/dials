#!/usr/bin/env python

from __future__ import division
from dials.scratch.dgw.source_model import source
from model_parameters import parameter, model_parameterisation
from scitbx import matrix
from math import pi
from dials.scratch.dgw.refinement \
    import dR_from_axis_and_angle, get_fd_gradients, random_param_shift

class source_parameterisation_orientation(model_parameterisation):
    '''implementation of parameterisation for the source (direction only)
    with angles expressed in mrad'''

    def __init__(self, source):

        # The state of the source model consists of the orientation of the
        # s0 vector that it is modelling. The initial state is a snapshot
        # of the s0 vector at the point of initialisation. Future states are
        # composed by rotations around axes perpendicular to that direction.
        #
        # The 'models' attribute refers to the source vector contained by this
        # model.

        ### Set up the initial state
        s0dir = matrix.col(source.get_s0().normalize().elems)
        istate = s0dir

        ### Set up the parameters
        s0_plane_dir1 = source.get_s0().ortho().normalize()
        s0_plane_dir2 = source.get_s0().cross(s0_plane_dir1).normalize()
        # rotation around s0_plane_dir1
        mu1 = parameter(.0, s0_plane_dir1, 'angle')
        # rotation around s0_plane_dir2
        mu2 = parameter(.0, s0_plane_dir2, 'angle')

        # build the parameter list in a specific,  maintained order
        p_list = [mu1, mu2]

        # set up the list of model objects being parameterised (here
        # just a single source model)
        models = [source]

        # set up the base class
        model_parameterisation.__init__(self, models, istate, p_list)

        # call compose to calculate all the derivatives
        self.compose()

    def compose(self):

        # extract direction from the initial state
        is0 = self._initial_state

        # extract parameters from the internal list
        mu1, mu2 = self._plist

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
        self._models[0].set_s0(s0mag * s0_new_dir)

        ### calculate derivatives of the state wrt parameters
        # derivative wrt mu1
        dMu21_dmu1 = Mu2 * dMu1_dmu1
        ds0_new_dir_dmu1 = dMu21_dmu1 * is0

        # derivative wrt mu2
        dMu21_dmu2 = dMu2_dmu2 * Mu1
        ds0_new_dir_dmu2 = dMu21_dmu2 * is0

        ### calculate derivatives of the attached source vector, converting
        ### parameters back to mrad, and store
        # derivative wrt mu1
        self._dstate_dp[0] = ds0_new_dir_dmu1 * s0mag / 1000.
        # derivative wrt mu2
        self._dstate_dp[1] = ds0_new_dir_dmu2 * s0mag / 1000.

        return

    def get_state(self):
        return matrix.col(self._models[0].get_s0())

if __name__ == '__main__':

    import random
    from libtbx.test_utils import approx_equal

    # make a random beam vector and parameterise it
    s0 = source(matrix.col.random(3, 0.5, 1.5))
    #s0p = source_parameterisation_orientation(s0)
    s0p = source_parameterisation_orientation(s0)

    # Let's do some basic tests. First, can we change parameter values and
    # update the modelled vector s0?
    s0_old = s0.get_s0()
    s0p.set_p([1000*0.1, 1000*0.1])
    assert(approx_equal(s0.get_s0().angle(s0_old), 0.1413033))

    # random initial orientations with a random parameter shift at each
    attempts = 1000
    failures = 0
    for i in range(attempts):

        # make a random beam vector and parameterise it
        s0 = source(matrix.col.random(3, 0.5, 1.5))
        s0p = source_parameterisation_orientation(s0)

        # apply a random parameter shift
        p_vals = s0p.get_p()
        p_vals = random_param_shift(p_vals, [1000*pi/9, 1000*pi/9])
        s0p.set_p(p_vals)

        # compare analytical and finite difference derivatives
        an_ds_dp = s0p.get_ds_dp()
        fd_ds_dp = get_fd_gradients(s0p, [1.e-5 * pi/180] * 2)

        for j in range(2):
            try:
                assert(approx_equal((fd_ds_dp[j] - an_ds_dp[j]),
                        matrix.col((0., 0., 0.)), eps = 1.e-6))
            except Exception:
                failures += 1
                print "for try", i
                print "failure for parameter number", j
                print "with fd_ds_dp = "
                print fd_ds_dp[j]
                print "and an_ds_dp = "
                print an_ds_dp[j]
                print "so that difference fd_ds_dp - an_ds_dp ="
                print fd_ds_dp[j] - an_ds_dp[j]

    if failures == 0: print "OK"
