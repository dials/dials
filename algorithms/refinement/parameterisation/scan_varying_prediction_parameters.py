#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# Python imports
from __future__ import division

# cctbx imports
from scitbx import matrix

# dials imports
from cctbx.array_family import flex
from dials_refinement_helpers_ext import *
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    DetectorSpacePredictionParameterisation

class VaryingCrystalPredictionParameterisation(DetectorSpacePredictionParameterisation):

    '''Support crystal parameterisations that vary with time (via its
    proxy of "observed image number"'''

    def _prepare(self):
        '''Cache required quantities that are not dependent on hkl'''

        # Same as _prepare for the parent class except we don't get
        # U and B from the model
        self._D = matrix.sqr(self._detector[0].get_D_matrix())
        self._s0 = matrix.col(self._beam.get_s0())
        self._axis = matrix.col(self._gonio.get_rotation_axis())

    def get_gradients(self, h, s, phi, obs_image_number):

        '''Adds obs_image_number for scan-varying parameters'''

        self._prepare()

        return self._get_gradients_core(h, s, phi, obs_image_number)

    def get_multi_gradients(self, match_list):
        '''
        Adds passing the observed image number to the gradient calc
        for scan-varying parameters
        '''

        self._prepare()

        # FIXME ObsPredMatch does not have an image attribute yet so
        # the following line will not work...
        return [self._get_gradients_core(m.H, m.Sc, m.Phic, m.image) for m in match_list]

    def _get_gradients_core(self, h, s, phi, obs_image_number):

        '''Calculate gradients of the prediction formula with respect to each
        of the parameters of the contained models, for reflection h with
        scattering vector s that reflects at rotation angle phi. That is,
        calculate dX/dp, dY/dp and dphi/dp. Scan-varying parameters (for
        the crystal) are evaluated at obs_image_number'''

        ### Calculate various quantities of interest for this reflection

        R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)

        # Get U and B for obs_image_number. Assume there is only one
        # parameterisation of each type
        xl_op = self._xl_orientation_parameterisations[0]
        xl_ucp = self._xl_unit_cell_parameterisations[0]
        U = xl_op.get_state(obs_image_number)
        B = xl_ucp.get_state(obs_image_number)

        # pv is the 'projection vector' for the reflection s.
        s = matrix.col(s)
        pv = self._D * s
        # r is the reciprocal lattice vector, in the lab frame
        r = R * U * B * h

        # All of the derivatives of phi have a common denominator, given by
        # (e X r).s0, where e is the rotation axis. Calculate this once, here.
        e_X_r = self._axis.cross(r)
        e_r_s0 = (e_X_r).dot(self._s0)

        try:
            assert abs(e_r_s0) > 1.e-6
        except AssertionError as e:
            print "(e X r).s0 too small:", e_r_s0
            print "for reflection", h
            print "with scattering vector", s
            print "where r =", r
            print "e =",matrix.col(self._gonio.get_rotation_axis())
            print "s0 =",self._s0
            print "U(t) =",U
            print "this reflection forms angle with the equatorial plane normal:"
            vecn = self._s0.cross(self._axis).normalize()
            print s.accute_angle(vecn)
            raise e
        # FIXME This is potentially dangerous! e_r_s0 -> 0 when the rotation
        # axis, beam vector and relp are coplanar. This occurs when a reflection
        # just touches the Ewald sphere. The derivatives of phi go to infinity
        # because any change moves it off this one position of grazing
        # incidence. How best to deal with this?

        ### Work through the parameterisations, calculating their contributions
        ### to derivatives d[pv]/dp and d[phi]/dp

        # Set up the lists of derivatives
        dpv_dp = []
        dphi_dp = []

        # Calculate derivatives of pv wrt each parameter of the FIRST detector
        # parameterisation only. All derivatives of phi are zero for detector
        # parameters
        if self._detector_parameterisations:
            self._detector_derivatives(dpv_dp, dphi_dp, pv)

        # Calc derivatives of pv and phi wrt each parameter of each beam
        # parameterisation that is present.
        if self._beam_parameterisations:
            self._beam_derivatives(dpv_dp, dphi_dp, r, e_X_r, e_r_s0)

        # Calc derivatives of pv and phi wrt each parameter of each
        # scan-varying crystal orientation parameterisation
        if self._xl_orientation_parameterisations:
            self._xl_orientation_derivatives(dpv_dp, dphi_dp, \
                    obs_image_number, B, R, h, s, e_X_r, e_r_s0)

        # Now derivatives of pv and phi wrt each parameter of each
        # scan-varying crystal unit cell parameterisation
        if self._xl_unit_cell_parameterisations:
            self._xl_unit_cell_derivatives(dpv_dp, dphi_dp, \
                    obs_image_number, U, R, h, s, e_X_r, e_r_s0)

        # calculate positional derivatives from d[pv]/dp
        pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e) for e in dpv_dp]
        dX_dp, dY_dp = zip(*pos_grad)

        return zip(dX_dp, dY_dp, dphi_dp)

    def _xl_orientation_derivatives(self, dpv_dp, dphi_dp, \
            obs_image_number, B, R, h, s, e_X_r, e_r_s0):

        '''Adds calculation at obs_image_number for scan-varying
        parameters'''

        for xlo in self._xl_orientation_parameterisations:
            dU_dxlo_p = xlo.get_ds_dp(obs_image_number)

            dr_dxlo_p = [R * dU_dxlo_p[i] * B * h for i in range(len(dU_dxlo_p))]

            dphi_dxlo_p = [- der.dot(s) / e_r_s0 for der in dr_dxlo_p]

            dpv_dxlo_p = [self._D * (dr_dxlo_p[i] + e_X_r * dphi_dxlo_p[i]) for i in range(len(dphi_dxlo_p))]

            dpv_dp.extend(dpv_dxlo_p)
            dphi_dp.extend(dphi_dxlo_p)

        return

    def _xl_unit_cell_derivatives(self, dpv_dp, dphi_dp, \
            obs_image_number, U, R, h, s, e_X_r, e_r_s0):

        '''Adds calculation at obs_image_number for scan-varying
        parameters'''

        for xluc in self._xl_unit_cell_parameterisations:
            dB_dxluc_p = xluc.get_ds_dp(obs_image_number)

            dr_dxluc_p = [R * U * dB_dxluc_p[i] * h for i
                              in range(len(dB_dxluc_p))]

            dphi_dxluc_p = [- der.dot(s) / e_r_s0 for der in dr_dxluc_p]

            dpv_dxluc_p = [self._D * (dr_dxluc_p[i] + e_X_r * dphi_dxluc_p[i]) for i in range(len(dr_dxluc_p))]

            dpv_dp.extend(dpv_dxluc_p)
            dphi_dp.extend(dphi_dxluc_p)

        return


    def _calc_dX_dp_and_dY_dp_from_dpv_dp(self, pv, der):
        '''helper function to calculate positional derivatives from dpv_dp using
        the quotient rule'''
        u = pv[0]
        v = pv[1]
        w = pv[2]
        w2 = w**2

        du_dp = der[0]
        dv_dp = der[1]
        dw_dp = der[2]

        dX_dp = du_dp / w - u * dw_dp / w2
        dY_dp = dv_dp / w - v * dw_dp / w2

        return dX_dp, dY_dp
