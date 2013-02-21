#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

from __future__ import division
from model_parameters import parameter, model_parameterisation
from dials.scratch.dgw.crystal_model import crystal
from scitbx import matrix
from dials.scratch.dgw.refinement \
    import get_fd_gradients, dR_from_axis_and_angle, random_param_shift
from math import pi
from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge

class crystal_orientation_parameterisation(model_parameterisation):
    '''A work-in-progress parameterisation for crystal orientation, with angles
    expressed in mrad'''

    def __init__(self, crystal):

        # The state of a crystal orientation parameterisation is an orientation
        # matrix '[U]'. The initial state is a snapshot of the crystal orientation
        # at the time of initialisation '[U0]'. Future states are composed by
        # rotations around axes of the phi-axis frame by Tait-Bryan angles.
        #
        # [U] = [Phi3][Phi2][Phi1][U0]

        ### Set up the initial state
        istate = crystal.get_U()

        ### Set up the parameters
        phi1 = parameter(.0, matrix.col((1, 0, 0)), 'angle')
        phi2 = parameter(.0, matrix.col((0, 1, 0)), 'angle')
        phi3 = parameter(.0, matrix.col((0, 0, 1)), 'angle')

        # build the parameter list in a specific,  maintained order
        p_list = [phi1, phi2, phi3]

        # set up the list of model objects being parameterised (here
        # just a single crystal model)
        models = [crystal]

        # set up the base class
        model_parameterisation.__init__(self, models, istate, p_list)

        # call compose to calculate all the derivatives
        self.compose()

    def compose(self):
        '''calculate state and derivatives'''

        # Extract orientation from the initial state
        U0 = self._initial_state

        # extract parameters from the internal list
        phi1, phi2, phi3 = self._plist

        # convert to radians
        phi1rad, phi2rad, phi3rad = (phi1.value / 1000., phi2.value / 1000.,
                                     phi3.value / 1000.)

        # compose rotation matrices and their first order derivatives
        Phi1 = (phi1.axis).axis_and_angle_as_r3_rotation_matrix(phi1rad, deg=False)
        dPhi1_dphi1 = dR_from_axis_and_angle(phi1.axis, phi1rad, deg=False)

        Phi2 = (phi2.axis).axis_and_angle_as_r3_rotation_matrix(phi2rad, deg=False)
        dPhi2_dphi2 = dR_from_axis_and_angle(phi2.axis, phi2rad, deg=False)

        Phi3 = (phi3.axis).axis_and_angle_as_r3_rotation_matrix(phi3rad, deg=False)
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

class crystal_unit_cell_parameterisation(model_parameterisation):
    '''A work-in-progress parameterisation for unit cell'''

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

        ### Set up the independent parameters
        p_list = [parameter(e) for e in X]

        # set up the list of model objects being parameterised (here
        # just a single crystal model)
        models = [crystal]

        # set up the base class
        model_parameterisation.__init__(self, models, istate, p_list)

        # call compose to calculate all the derivatives
        self.compose()

    def compose(self):
        '''calculate state and derivatives'''

        p_vals = [p.value for p in self._plist]

        # set parameter values in the symmetrizing object and obtain new B
        newB = matrix.sqr(self._S.backward_orientation(p_vals).reciprocal_matrix())

        # Now pass new B to the crystal model
        self._models[0].set_B(newB)

        # returns the independent parameters given the set_orientation() B matrix
        # used here for side effects
        self._S.forward_independent_parameters()

        # get the gradients
        self._dstate_dp = [matrix.sqr(e) for e in self._S.forward_gradients()]

        return

    def get_state(self):
        return matrix.sqr(self._models[0].get_B())

if __name__ == '__main__':

    import random
    from libtbx.test_utils import approx_equal
    from cctbx.uctbx import unit_cell

    def random_direction_close_to(vector):
        return vector.rotate_around_origin(matrix.col(
                    (random.random(),
                     random.random(),
                     random.random())).normalize(),
                     random.gauss(0, 1.0),  deg = True)

    # make a random crystal and parameterise it
    a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
    b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
    c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
    xl = crystal(a, b, c)

    xl_op = crystal_orientation_parameterisation(xl)
    xl_ucp = crystal_unit_cell_parameterisation(xl)

    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))

    # compare analytical and finite difference derivatives
    an_ds_dp = xl_op.get_ds_dp()
    fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * 3)
    for e, f in zip(an_ds_dp, fd_ds_dp):
        assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    an_ds_dp = xl_ucp.get_ds_dp()
    fd_ds_dp = get_fd_gradients(xl_ucp, [1.e-7] * xl_ucp.num_free())
    for e, f in zip(an_ds_dp, fd_ds_dp):
        print e
        print f
        assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    # random initial orientations with a random parameter shift at each
    attempts = 100
    failures = 0
    for i in range(attempts):

        # make a random crystal and parameterise it
        a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
        b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
        c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))
        xl = crystal(a, b, c)
        xl_op = crystal_orientation_parameterisation(xl)
        xl_uc = crystal_unit_cell_parameterisation(xl)

        # apply a random parameter shift to the orientation
        p_vals = xl_op.get_p()
        p_vals = random_param_shift(p_vals, [1000*pi/9, 1000*pi/9,
                                             1000*pi/9])
        xl_op.set_p(p_vals)

        # compare analytical and finite difference derivatives
        xl_op_an_ds_dp = xl_op.get_ds_dp()
        xl_op_fd_ds_dp = get_fd_gradients(xl_op, [1.e-5 * pi/180] * 3)

        # apply a random parameter shift to the unit cell

        print "\nCYCLE", i, "\n"
        print "apply random parameter shift"
        p_vals = xl_uc.get_p()
        cell_params = xl.get_unit_cell().parameters()
        print "old unit cell",cell_params
        cell_params = random_param_shift(cell_params, [1.] * 6)
        new_uc = unit_cell(cell_params)
        print "new unit cell",cell_params
        newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
        S = symmetrize_reduce_enlarge(xl.get_space_group())
        S.set_orientation(orientation=newB)
        X = S.forward_independent_parameters()
        print "old p_vals", p_vals
        print "new p_vals", X, "\n"
        print "set these parameters in the model parameterisation"
        xl_uc.set_p(X)

        xl_uc_an_ds_dp = xl_ucp.get_ds_dp()
        print "\nnow doing finite differences about each parameter in turn"
        xl_uc_fd_ds_dp = get_fd_gradients(xl_ucp, [1.e-7] * xl_ucp.num_free())

        for j in range(3):
            try:
                assert(approx_equal((xl_op_fd_ds_dp[j] - xl_op_an_ds_dp[j]),
                                    null_mat, eps = 1.e-6))
            except Exception:
                failures += 1
                print "for try", i
                print "failure for parameter number", j
                print "of the orientation parameterisation"
                print "with fd_ds_dp = "
                print fd_ds_dp[j]
                print "and an_ds_dp = "
                print an_ds_dp[j]
                print "so that difference fd_ds_dp - an_ds_dp ="
                print fd_ds_dp[j] - an_ds_dp[j]

        for j in range(xl_ucp.num_free()):
            try:
                assert(approx_equal((xl_uc_fd_ds_dp[j] - xl_uc_an_ds_dp[j]),
                                    null_mat, eps = 1.e-6))
            except Exception:
                failures += 1
                print "for try", i
                print "failure for parameter number", j
                print "of the unit cell parameterisation"
                print "with fd_ds_dp = "
                print fd_ds_dp[j]
                print "and an_ds_dp = "
                print an_ds_dp[j]
                print "so that difference fd_ds_dp - an_ds_dp ="
                print fd_ds_dp[j] - an_ds_dp[j]

    if failures == 0: print "OK"
