#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

#### Python and general cctbx imports

from __future__ import division
from scitbx import matrix

#### Import model parameterisations

from dials.array_family import flex
from dials_refinement_helpers_ext import *

from dials.algorithms.refinement.parameterisation import \
    PredictionParameterisation

class PredictionParameterisationXY(PredictionParameterisation):
    '''
    Differs from PredictionParameterisation only by excluding the
    goniometer model
    '''

    def __init__(self,
                 detector_model,
                 beam_model,
                 crystal_model,
                 detector_parameterisations = None,
                 beam_parameterisations = None,
                 xl_orientation_parameterisations = None,
                 xl_unit_cell_parameterisations = None):

        # References to the underlying models
        self._detector = detector_model
        self._beam = beam_model
        self._crystal = crystal_model

        # Keep references to all parameterised models
        self._detector_parameterisations = detector_parameterisations
        self._beam_parameterisations = beam_parameterisations
        self._xl_orientation_parameterisations = \
            xl_orientation_parameterisations
        self._xl_unit_cell_parameterisations = \
            xl_unit_cell_parameterisations

        self._length = self._len()

        # Fill out remaining attributes by a call to prepare
        self.prepare()

    def prepare(self):
        '''Cache required quantities that are not dependent on hkl'''

        # Obtain various quantities of interest from the experimental model
        self._D_mats = [matrix.sqr(p.get_D_matrix()) for p in self._detector]
        self._s0 = matrix.col(self._beam.get_s0())
        self._U = self._crystal.get_U()
        self._B = self._crystal.get_B()
        self._UB = self._U * self._B
        #self._axis = matrix.col(self._gonio.get_rotation_axis())

    def get_gradients(self, h, s, phi, panel_id, obs_image_number = None):
        '''
        Calculate gradients of the prediction formula with respect to each
        of the parameters of the contained models, for the reflection with
        scattering vector s that intersects panel with panel_id.

        To be implemented by a derived class, which determines the space of the
        prediction formula (e.g. we calculate dX/dp, dY/dp, dphi/dp for the
        prediction formula expressed in detector space, but components of
        d\vec{r}/dp for the prediction formula in reciprocal space

        obs_image_number included to match the interface of a scan-
        varying version of the class
        '''

        # extract the right panel matrix
        self._D = self._D_mats[panel_id]

        return self._get_gradients_core(h, s, phi, panel_id)


class DetectorSpacePredictionParameterisation(PredictionParameterisationXY):
    '''
    Concrete class that inherits functionality of the
    PredictionParameterisation parent class and provides a detector space
    implementation of the get_gradients function.

    Untested for multiple sensor detectors.
    '''

    def _get_gradients_core(self, h, s, phi, panel_id):

        '''Calculate gradients of the prediction formula with respect to
        each of the parameters of the contained models, for reflection h
        that reflects at rotation angle phi with scattering vector s that
        intersects panel panel_id. That is, calculate dX/dp, dY/dp and
        dphi/dp'''

        ### Calculate various quantities of interest for this reflection

        R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)

        # pv is the 'projection vector' for the reflection s.
        s = matrix.col(s)
        pv = self._D * s
        # r is the reciprocal lattice vector, in the lab frame
        r = R * self._UB * h

        # All of the derivatives of phi have a common denominator, given by
        # (e X r).s0, where e is the rotation axis. Calculate this once, here.
        #e_X_r = self._axis.cross(r)
        #e_r_s0 = (e_X_r).dot(self._s0)

        # Note relationship between e_r_s0 and zeta_factor. Uncommenting the
        # code below shows that s0.(e X r) = zeta * |s X s0|
        #from dials.algorithms.reflection_basis import zeta_factor
        #from libtbx.test_utils import approx_equal
        #z = zeta_factor(self._axis, self._s0, s)
        #ss0 = (s.cross(self._s0)).length()
        #assert approx_equal(e_r_s0, z * ss0)

        #try:
        #    assert abs(e_r_s0) > 1.e-6
        #except AssertionError as e:
        #    print "(e X r).s0 too small:", e_r_s0
        #    print "for reflection", h
        #    print "with scattering vector", s
        #    print "where r =", r
        #    print "e =",matrix.col(self._gonio.get_rotation_axis())
        #    print "s0 =",self._s0
        #    print "U =",self._U
        #    print "this reflection forms angle with the equatorial plane normal:"
        #    vecn = self._s0.cross(self._axis).normalize()
        #    print s.accute_angle(vecn)
        #    raise e
        # This is potentially dangerous! e_r_s0 -> 0 when the rotation
        # axis, beam vector and relp are coplanar. This occurs when a reflection
        # just touches the Ewald sphere.

        ### Work through the parameterisations, calculating their contributions
        ### to derivatives d[pv]/dp and d[phi]/dp

        # Set up the lists of derivatives
        dpv_dp = []
        #dphi_dp = []

        # Calculate derivatives of pv wrt each parameter of the FIRST detector
        # parameterisation only. All derivatives of phi are zero for detector
        # parameters
        if self._detector_parameterisations:
            self._detector_derivatives(dpv_dp, pv, panel_id)

        # Calc derivatives of pv and phi wrt each parameter of each beam
        # parameterisation that is present.
        if self._beam_parameterisations:
            self._beam_derivatives(dpv_dp, r)

        # Calc derivatives of pv and phi wrt each parameter of each crystal
        # orientation parameterisation that is present.
        if self._xl_orientation_parameterisations:
            self._xl_orientation_derivatives(dpv_dp, R, h)

        # Now derivatives of pv and phi wrt each parameter of each crystal unit
        # cell parameterisation that is present.
        if self._xl_unit_cell_parameterisations:
            self._xl_unit_cell_derivatives(dpv_dp, R, h)

        # calculate positional derivatives from d[pv]/dp
        pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e) for e in dpv_dp]
        dX_dp, dY_dp = zip(*pos_grad)

        return zip(dX_dp, dY_dp)

    def _detector_derivatives(self, dpv_dp, pv, panel_id):

        '''helper function to extend the derivatives lists by
        derivatives of the detector parameterisations'''

        for idet, det in enumerate(self._detector_parameterisations):
            if idet == 0:
                dd_ddet_p = det.get_ds_dp(multi_state_elt=panel_id)
                dpv_ddet_p = [- self._D * e * pv for e in dd_ddet_p]
            else:
                dpv_ddet_p = [matrix.col((0., 0., 0.))] * len(dd_ddet_p)

            #dphi_ddet_p = [0.] * len(dd_ddet_p)

            dpv_dp.extend(dpv_ddet_p)
            #dphi_dp.extend(dphi_ddet_p)

        return

    def _beam_derivatives(self, dpv_dp, r):

        '''helper function to extend the derivatives lists by
        derivatives of the beam parameterisations'''

        for src in self._beam_parameterisations:
            ds0_dsrc_p = src.get_ds_dp()
            #dphi_dsrc_p = [- r.dot(ds0_dsrc_p[i]) / e_r_s0 for i
            #                  in range(len(ds0_dsrc_p))]
            dpv_dsrc_p = [self._D * e for e in ds0_dsrc_p]

            dpv_dp.extend(dpv_dsrc_p)
            #dphi_dp.extend(dphi_dsrc_p)

            return

    def _xl_orientation_derivatives(self, dpv_dp, R, h):

        '''helper function to extend the derivatives lists by
        derivatives of the crystal orientation parameterisations'''

        for xlo in self._xl_orientation_parameterisations:
            dU_dxlo_p = xlo.get_ds_dp()

            dr_dxlo_p = [R * e * self._B * h for e in dU_dxlo_p]

            #dphi_dxlo_p = [- der.dot(s) / e_r_s0 for der in dr_dxlo_p]

            dpv_dxlo_p = [self._D * e for e in dr_dxlo_p]

            dpv_dp.extend(dpv_dxlo_p)
            #dphi_dp.extend(dphi_dxlo_p)

        return

    def _xl_orientation_derivatives(self, dpv_dp, R, h):

        '''helper function to extend the derivatives lists by
        derivatives of the crystal orientation parameterisations'''

        for xlo in self._xl_orientation_parameterisations:
            dU_dxlo_p = xlo.get_ds_dp()

            dr_dxlo_p = [R * e * self._B * h for e in dU_dxlo_p]

            #dphi_dxlo_p = [- der.dot(s) / e_r_s0 for der in dr_dxlo_p]

            dpv_dxlo_p = [self._D * e for e in dr_dxlo_p]

            dpv_dp.extend(dpv_dxlo_p)
            #dphi_dp.extend(dphi_dxlo_p)

        return

    def _xl_unit_cell_derivatives(self, dpv_dp, R, h):

        '''helper function to extend the derivatives lists by
        derivatives of the crystal unit cell parameterisations'''

        for xluc in self._xl_unit_cell_parameterisations:
            dB_dxluc_p = xluc.get_ds_dp()

            dr_dxluc_p = [R * self._U * e * h for e in dB_dxluc_p]

            #phi_dxluc_p = [- der.dot(s) / e_r_s0 for der in dr_dxluc_p]

            dpv_dxluc_p = [self._D * e for e in dr_dxluc_p]

            dpv_dp.extend(dpv_dxluc_p)
            #dphi_dp.extend(dphi_dxluc_p)

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
