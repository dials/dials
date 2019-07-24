#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    PredictionParameterisation,
    SparseGradientVectorMixin,
)
from dials_refinement_helpers_ext import dRq_de


class StillsPredictionParameterisation(PredictionParameterisation):
    """Concrete class that inherits functionality of the
    PredictionParameterisation parent class and provides a detector space
    implementation of the get_gradients function for still images. Gradients of
    the minimum rotation to the Ewald sphere, DeltaPsi, are calculated for use
    as a restraint in refinement"""

    _grad_names = ("dX_dp", "dY_dp", "dDeltaPsi_dp")

    def __init__(self, *args, **kwargs):
        super(StillsPredictionParameterisation, self).__init__(*args, **kwargs)
        # check that a goniometer parameterisation is not passed in
        assert not self._goniometer_parameterisations
        return

    def _local_setup(self, reflections):
        """Setup additional attributes used in gradients calculation. These are
        specific to the stills-type prediction parameterisations that rotates the
        relp q onto the Ewald sphere by an angle DeltaPsi before prediction"""

        self._DeltaPsi = reflections["delpsical.rad"]

        # q is the reciprocal lattice vector, in the lab frame
        self._q = self._UB * self._h
        self._q_scalar = self._q.norms()
        self._qq = self._q_scalar * self._q_scalar

        # r is the reciprocal lattice vector rotated to the Ewald sphere
        self._r = self._s1 - self._s0

        # we also need the unit directions q0 and s0u
        self._q0 = self._q.each_normalize()
        self._s0u = self._s0.each_normalize()

        # e1 is the unit vector about which DeltaPsi rotation is defined
        self._e1 = self._q0.cross(self._s0u).each_normalize()

        # q1 completes an orthonormal set with q0 and e1
        self._q1 = self._q0.cross(self._e1).each_normalize()

        # c0 completes an orthonormal set with s0u and e1
        self._c0 = self._s0u.cross(self._e1).each_normalize()

        # we want the wavelength
        self._wavelength = 1.0 / self._s0.norms()

        return

    def _beam_derivatives(self, isel, parameterisation=None, reflections=None):
        """helper function to extend the derivatives lists by derivatives of the
        beam parameterisations"""

        # Get required data
        s0 = self._s0.select(isel)
        s0u = self._s0u.select(isel)
        wl = self._wavelength.select(isel)
        r = self._r.select(isel)
        e1 = self._e1.select(isel)
        q = self._q.select(isel)
        c0 = self._c0.select(isel)
        DeltaPsi = self._DeltaPsi.select(isel)
        D = self._D.select(isel)

        # get the derivatives of the beam vector wrt the parameters
        ds0_dbeam_p = parameterisation.get_ds_dp(use_none_as_null=True)

        dDeltaPsi_dp = []
        dpv_dp = []

        # loop through the parameters
        for der in ds0_dbeam_p:

            if der is None:
                dpv_dp.append(None)
                dDeltaPsi_dp.append(None)
                continue

            # repeat the derivative in an array
            ds0 = flex.vec3_double(len(s0u), der.elems)

            # we need the derivative of the unit beam direction too. This requires
            # scaling by the wavelength and projection onto the Ewald sphere
            scaled = ds0 * wl
            ds0u = scaled.dot(c0) * c0 + scaled.dot(e1) * e1

            # calculate the derivative of DeltaPsi for this parameter
            dDeltaPsi = -1.0 * (r.dot(ds0)) / (e1.cross(r).dot(s0))
            dDeltaPsi_dp.append(dDeltaPsi)

            # calculate (d[r]/d[e1])(d[e1]/dp)
            de1_dp = c0.cross(ds0u)
            dr_de1 = dRq_de(DeltaPsi, e1, q)
            drde_dedp = dr_de1 * de1_dp

            # dp = 1.e-8 # finite step size for the parameter
            # del_e1 = de1_dp * dp
            # e1f = e1 + del_e1 * 0.5
            # rfwd = q.rotate_around_origin(e1f, DeltaPsi)
            # e1r = e1 - del_e1 * 0.5
            # rrev = q.rotate_around_origin(e1r, DeltaPsi)
            # drde_dedp = (rfwd - rrev) * (1 / dp)

            # calculate the derivative of pv for this parameter
            dpv_dp.append(D * (ds0 + e1.cross(r) * dDeltaPsi + drde_dedp))

        return dpv_dp, dDeltaPsi_dp

    def _xl_orientation_derivatives(
        self, isel, parameterisation=None, reflections=None
    ):
        """helper function to extend the derivatives lists by derivatives of the
        crystal orientation parameterisations"""

        # Get required data
        B = self._B.select(isel)
        h = self._h.select(isel)
        e1 = self._e1.select(isel)
        DeltaPsi = self._DeltaPsi.select(isel)
        s1 = self._s1.select(isel)
        q = self._q.select(isel)
        q_scalar = self._q_scalar.select(isel)
        qq = self._qq.select(isel)
        q0 = self._q0.select(isel)
        r = self._r.select(isel)
        s0 = self._s0.select(isel)
        s0u = self._s0u.select(isel)
        D = self._D.select(isel)

        # get derivatives of the U matrix wrt the parameters
        dU_dxlo_p = parameterisation.get_ds_dp(use_none_as_null=True)

        dDeltaPsi_dp = []
        dpv_dp = []

        # loop through the parameters
        for der in dU_dxlo_p:

            if der is None:
                dpv_dp.append(None)
                dDeltaPsi_dp.append(None)
                continue

            der_mat = flex.mat3_double(len(B), der.elems)

            # calculate the derivative of q for this parameter
            dq = der_mat * B * h

            # calculate the derivative of r for this parameter
            dr = dq.rotate_around_origin(e1, DeltaPsi)

            # calculate the derivative of DeltaPsi for this parameter
            dDeltaPsi = -1.0 * (dr.dot(s1)) / (e1.cross(r).dot(s0))
            dDeltaPsi_dp.append(dDeltaPsi)

            # derivative of the axis e1
            q_dot_dq = q.dot(dq)
            dq0 = (q_scalar * dq - (q_dot_dq * q0)) / qq
            de1_dp = dq0.cross(s0u)

            # calculate (d[r]/d[e1])(d[e1]/dp)
            # dr_de1 = dRq_de(DeltaPsi, e1, q)
            # drde_dedp = dr_de1 * de1_dp

            # The above calculation is not correct. Instead, calculate the partial
            # derivative of r wrt change in rotation axis e1 by finite differences
            dp = 1.0e-8
            del_e1 = de1_dp * dp
            e1f = e1 + del_e1 * 0.5
            rfwd = q.rotate_around_origin(e1f, DeltaPsi)
            e1r = e1 - del_e1 * 0.5
            rrev = q.rotate_around_origin(e1r, DeltaPsi)
            drde_dedp = (rfwd - rrev) * (1 / dp)

            # calculate the derivative of pv for this parameter
            dpv_dp.append(D * (dr + e1.cross(r) * dDeltaPsi + drde_dedp))

        return dpv_dp, dDeltaPsi_dp

    def _xl_unit_cell_derivatives(self, isel, parameterisation=None, reflections=None):
        """helper function to extend the derivatives lists by derivatives of the
        crystal unit cell parameterisations"""

        # Get required data
        U = self._U.select(isel)
        h = self._h.select(isel)
        e1 = self._e1.select(isel)
        DeltaPsi = self._DeltaPsi.select(isel)
        s1 = self._s1.select(isel)
        q = self._q.select(isel)
        q_scalar = self._q_scalar.select(isel)
        qq = self._qq.select(isel)
        q0 = self._q0.select(isel)
        r = self._r.select(isel)
        s0 = self._s0.select(isel)
        s0u = self._s0u.select(isel)
        D = self._D.select(isel)

        # get derivatives of the B matrix wrt the parameters
        dB_dxluc_p = parameterisation.get_ds_dp(use_none_as_null=True)

        dDeltaPsi_dp = []
        dpv_dp = []

        # loop through the parameters
        for der in dB_dxluc_p:

            if der is None:
                dpv_dp.append(None)
                dDeltaPsi_dp.append(None)
                continue

            der_mat = flex.mat3_double(len(U), der.elems)

            # calculate the derivative of q for this parameter
            dq = U * der_mat * h

            # calculate the derivative of r for this parameter
            dr = dq.rotate_around_origin(e1, DeltaPsi)

            # calculate the derivative of DeltaPsi for this parameter
            dDeltaPsi = -1.0 * (dr.dot(s1)) / (e1.cross(r).dot(s0))
            dDeltaPsi_dp.append(dDeltaPsi)

            # derivative of the axis e1
            q_dot_dq = q.dot(dq)
            dq0 = (q_scalar * dq - (q_dot_dq * q0)) / qq
            de1_dp = dq0.cross(s0u)

            # calculate (d[r]/d[e1])(d[e1]/dp)
            # dr_de1 = dRq_de(DeltaPsi, e1, q)
            # drde_dedp = dr_de1 * de1_dp

            # The above calculation is not correct. Instead, calculate the partial
            # derivative of r wrt change in rotation axis e1 by finite differences
            dp = 1.0e-8
            del_e1 = de1_dp * dp
            e1f = e1 + del_e1 * 0.5
            rfwd = q.rotate_around_origin(e1f, DeltaPsi)
            e1r = e1 - del_e1 * 0.5
            rrev = q.rotate_around_origin(e1r, DeltaPsi)
            drde_dedp = (rfwd - rrev) * (1 / dp)

            # calculate the derivative of pv for this parameter
            dpv = D * (dr + e1.cross(r) * dDeltaPsi + drde_dedp)
            dpv_dp.append(dpv)

        return dpv_dp, dDeltaPsi_dp

    @staticmethod
    def _calc_dX_dp_and_dY_dp_from_dpv_dp(w_inv, u_w_inv, v_w_inv, dpv_dp):
        """helper function to calculate positional derivatives from
        dpv_dp using the quotient rule"""

        dX_dp = []
        dY_dp = []

        for der in dpv_dp:
            if der is None:
                dX_dp.append(None)
                dY_dp.append(None)
            else:
                du_dp, dv_dp, dw_dp = der.parts()

                dX_dp.append(w_inv * (du_dp - dw_dp * u_w_inv))
                dY_dp.append(w_inv * (dv_dp - dw_dp * v_w_inv))

        return dX_dp, dY_dp


class StillsPredictionParameterisationSparse(
    SparseGradientVectorMixin, StillsPredictionParameterisation
):
    """A version of StillsPredictionParameterisation that uses a sparse matrix
    data structure for memory efficiency when there are a large number of
    Experiments"""

    pass


class SphericalRelpStillsPredictionParameterisation(StillsPredictionParameterisation):
    """Modified StillsPredictionParameterisation for the model that assumes
    relps are spherical and prediction requires that this sphere intersects the
    Ewald sphere, not that the relp centre is rotated onto the Ewald sphere.
    Gradients of DeltaPsi are still calculated for use as a restraint."""

    def _local_setup(self, reflections):
        """Setup additional attributes used in gradients calculation. These are
        specific to the stills-type prediction parameterisations that leaves the
        relp q generally off the Ewald sphere but assumes it has spherical extent
        which intersects the Ewald sphere. Gradients of the minimum rotation
        angle DeltaPsi are still calculated for use as a restraint"""

        self._DeltaPsi = reflections["delpsical.rad"]

        # q is the reciprocal lattice vector, in the lab frame
        self._q = self._UB * self._h

        # quantities involving the direct beam vector
        self._q_s0 = self._q + self._s0
        s = self._q_s0.norms()
        ss = self._q_s0.dot(self._q_s0)
        sss = s * ss
        self._inv_s = 1.0 / s
        self._inv_sss = 1.0 / sss

        # r is the reciprocal lattice vector rotated to the Ewald sphere, required
        # for derivatives of DeltaPsi
        self._r = self._s1 - self._s0

        # we also need the unit directions q0 and s0u
        self._q0 = self._q.each_normalize()
        self._s0u = self._s0.each_normalize()

        # e1 is the unit vector about which DeltaPsi rotation is defined
        self._e1 = self._q0.cross(self._s0u).each_normalize()

        # q1 completes an orthonormal set with q0 and e1
        self._q1 = self._q0.cross(self._e1).each_normalize()

        # c0 completes an orthonormal set with s0u and e1
        self._c0 = self._s0u.cross(self._e1).each_normalize()

        # we want the wavenumber and wavelength
        self._nu = self._s0.norms()
        self._wavelength = 1.0 / self._nu

        return

    def _beam_derivatives(self, isel, parameterisation=None, reflections=None):
        """helper function to extend the derivatives lists by derivatives of the
        beam parameterisations"""

        # Get required data
        s0 = self._s0.select(isel)
        s0u = self._s0u.select(isel)
        nu = self._nu.select(isel)
        r = self._r.select(isel)
        e1 = self._e1.select(isel)
        q_s0 = self._q_s0.select(isel)
        inv_s = self._inv_s.select(isel)
        inv_sss = self._inv_sss.select(isel)
        D = self._D.select(isel)

        # get the derivatives of the beam vector wrt the parameters
        ds0_dbeam_p = parameterisation.get_ds_dp(use_none_as_null=True)

        dDeltaPsi_dp = []
        dpv_dp = []

        # loop through the parameters
        for der in ds0_dbeam_p:

            if der is None:
                dpv_dp.append(None)
                dDeltaPsi_dp.append(None)
                continue

            # repeat the derivative in an array
            ds0 = flex.vec3_double(len(s0u), der.elems)

            # we need the derivative of the unit beam direction too. This requires
            # scaling by the wavelength and projection onto the Ewald sphere
            # scaled = ds0 * wl
            # ds0u = scaled.dot(c0) * c0 + scaled.dot(e1) * e1

            # calculate the derivative of DeltaPsi for this parameter
            dDeltaPsi = -1.0 * (r.dot(ds0)) / (e1.cross(r).dot(s0))
            dDeltaPsi_dp.append(dDeltaPsi)

            # term 1
            term1 = (s0u.dot(ds0) * q_s0 + nu * ds0) * inv_s

            # term 2
            term2 = (nu * q_s0 * q_s0.dot(ds0)) * inv_sss

            # calculate the derivative of pv for this parameter
            dpv_dp.append(D * (term1 - term2))

        return dpv_dp, dDeltaPsi_dp

    def _xl_orientation_derivatives(
        self, isel, parameterisation=None, reflections=None
    ):
        """helper function to extend the derivatives lists by
        derivatives of the crystal orientation parameterisations"""

        # Get required data
        B = self._B.select(isel)
        h = self._h.select(isel)
        e1 = self._e1.select(isel)
        DeltaPsi = self._DeltaPsi.select(isel)
        s1 = self._s1.select(isel)
        nu = self._nu.select(isel)
        q_s0 = self._q_s0.select(isel)
        inv_s = self._inv_s.select(isel)
        inv_sss = self._inv_sss.select(isel)
        r = self._r.select(isel)
        s0 = self._s0.select(isel)
        D = self._D.select(isel)

        # get derivatives of the U matrix wrt the parameters
        dU_dxlo_p = parameterisation.get_ds_dp(use_none_as_null=True)

        dDeltaPsi_dp = []
        dpv_dp = []

        # loop through the parameters
        for der in dU_dxlo_p:

            if der is None:
                dpv_dp.append(None)
                dDeltaPsi_dp.append(None)
                continue

            der_mat = flex.mat3_double(len(B), der.elems)

            # calculate the derivative of q for this parameter
            dq = der_mat * B * h

            # calculate the derivative of r for this parameter
            dr = dq.rotate_around_origin(e1, DeltaPsi)

            # calculate the derivative of DeltaPsi for this parameter
            dDeltaPsi = -1.0 * (dr.dot(s1)) / (e1.cross(r).dot(s0))
            dDeltaPsi_dp.append(dDeltaPsi)

            # term 1
            term1 = (nu * dq) * inv_s

            # term 2
            term2 = (nu * q_s0 * q_s0.dot(dq)) * inv_sss

            # calculate the derivative of pv for this parameter
            dpv_dp.append(D * (term1 - term2))

        return dpv_dp, dDeltaPsi_dp

    def _xl_unit_cell_derivatives(self, isel, parameterisation=None, reflections=None):
        """helper function to extend the derivatives lists by
        derivatives of the crystal unit cell parameterisations"""

        # Get required data
        U = self._U.select(isel)
        h = self._h.select(isel)
        e1 = self._e1.select(isel)
        DeltaPsi = self._DeltaPsi.select(isel)
        s1 = self._s1.select(isel)
        nu = self._nu.select(isel)
        q_s0 = self._q_s0.select(isel)
        inv_s = self._inv_s.select(isel)
        inv_sss = self._inv_sss.select(isel)
        r = self._r.select(isel)
        s0 = self._s0.select(isel)
        D = self._D.select(isel)

        # get derivatives of the B matrix wrt the parameters
        dB_dxluc_p = parameterisation.get_ds_dp(use_none_as_null=True)

        dDeltaPsi_dp = []
        dpv_dp = []

        # loop through the parameters
        for der in dB_dxluc_p:

            if der is None:
                dpv_dp.append(None)
                dDeltaPsi_dp.append(None)
                continue

            der_mat = flex.mat3_double(len(U), der.elems)

            # calculate the derivative of q for this parameter
            dq = U * der_mat * h

            # calculate the derivative of r for this parameter
            dr = dq.rotate_around_origin(e1, DeltaPsi)

            # calculate the derivative of DeltaPsi for this parameter
            dDeltaPsi = -1.0 * (dr.dot(s1)) / (e1.cross(r).dot(s0))
            dDeltaPsi_dp.append(dDeltaPsi)

            # term 1
            term1 = (nu * dq) * inv_s

            # term 2
            term2 = (nu * q_s0 * q_s0.dot(dq)) * inv_sss

            # calculate the derivative of pv for this parameter
            dpv = D * (term1 - term2)
            dpv_dp.append(dpv)

        return dpv_dp, dDeltaPsi_dp


class SphericalRelpStillsPredictionParameterisationSparse(
    SparseGradientVectorMixin, SphericalRelpStillsPredictionParameterisation
):
    """A version of SphericalRelpStillsPredictionParameterisation that uses a
    sparse matrix data structure for memory efficiency when there are a large
    number of Experiments"""

    pass
