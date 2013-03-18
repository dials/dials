# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

#### Python and general cctbx imports

from __future__ import division
from scitbx import matrix

#### Import model parameterisations

from dials.scratch.dgw.refinement.detector_parameters import \
    detector_parameterisation_single_sensor
from dials.scratch.dgw.refinement.source_parameters import \
    beam_parameterisation_orientation
from dials.scratch.dgw.refinement.crystal_parameters import \
    crystal_orientation_parameterisation, crystal_unit_cell_parameterisation
from cctbx.array_family import flex
from dials_refinement_ext import *

class prediction_parameterisation(object):
    '''
    Abstract interface for a class that groups together model parameterisations
    relating to diffraction geometry and provides:

    * A list of all free parameters concatenated from each of the models, with a
      getter and setter method that delegates to the contained models
    * Derivatives of the reflection prediction equation with respect to each of
      these free parameters

    Derived classes determine whether the reflection prediction equation is
    expressed in detector space (X, Y, phi) or orthogonalised reciprocal space.

    It is assumed that the provided model parameterisations will be one of four
    types:

    * Detector parameterisations
    * Beam parameterisations
    * Crystal orientation parameterisations
    * Crystal unit cell parameterisations

    Not all of these must be supplied, and there may be multiple instances
    of one type (independently parameterised detectors, for example). This
    implies the type of the constructor arguments are either list (of one or
    more elements) or None.

    We also need access to the underlying models that are parameterised. The
    model parameterisation objects do not provide access to these models,
    as that is not their job. Instead we construct this object with direct
    access to each of the models. At this stage we assume only one of each model
    can exist, and *must* be present (otherwise reflection prediction does not
    make sense). So, we need each of:

    * A detector model (single point of reference for every sensor in the
      experiment)
    * A beam model
    * A crystal model
    * A goniometer model (not yet parameterised, but required for the equations)

    A class implementing prediction_parameterisation is used by a refinery
    object directly, which takes the list of parameters, and indirectly via a
    target function object, which takes the list of derivatives and composes the
    derivatives of a target function from them.
    '''

    def __init__(self,
                 detector_model,
                 beam_model,
                 crystal_model,
                 goniometer_model,
                 detector_parameterisations = None,
                 beam_parameterisations = None,
                 crystal_orientation_parameterisations = None,
                 crystal_unit_cell_parameterisations = None):

        # References to the underlying models
        self._detector = detector_model
        self._beam = beam_model
        self._crystal = crystal_model
        self._gonio = goniometer_model

        # Sanity checks
        #if detector_parameterisations:
        #    for model in detector_parameterisations:
        #        # TODO replace detector_parameterisation_single_sensor with a
        #        # general multi sensor detector_parameterisation when available
        #        assert isinstance(
        #            model, detector_parameterisation_single_sensor)
        #
        #if beam_parameterisations:
        #    for model in beam_parameterisations:
        #        assert isinstance(model, beam_parameterisation_orientation)
        #
        #if crystal_orientation_parameterisations:
        #    for model in crystal_orientation_parameterisations:
        #        assert isinstance(model, crystal_orientation_parameterisation)
        #
        #if crystal_unit_cell_parameterisations:
        #    for model in crystal_unit_cell_parameterisations:
        #        assert isinstance(model, crystal_unit_cell_parameterisation)

        # Keep references to all parameterised models
        self._detector_parameterisations = detector_parameterisations
        self._beam_parameterisations = beam_parameterisations
        self._xl_orientation_parameterisations = \
            crystal_orientation_parameterisations
        self._xl_unit_cell_parameterisations = \
            crystal_unit_cell_parameterisations

        self._length = self._len()

    def _len(self):
        length = 0
        if self._detector_parameterisations:
            for model in self._detector_parameterisations:
                length += model.num_free()

        if self._beam_parameterisations:
            for model in self._beam_parameterisations:
                length += model.num_free()

        if self._xl_orientation_parameterisations:
            for model in self._xl_orientation_parameterisations:
                length += model.num_free()

        if self._xl_unit_cell_parameterisations:
            for model in self._xl_unit_cell_parameterisations:
                length += model.num_free()

        return length

    def __len__(self):
        return self._length

    def get_p(self):
        '''return a concatenated list of parameters from each of the components
        in the global model'''

        global_p_list = []
        if self._detector_parameterisations:
            det_plists = [x.get_p() for x in self._detector_parameterisations]
            params = [x for l in det_plists for x in l]
            global_p_list.extend(params)

        if self._beam_parameterisations:
            src_plists = [x.get_p() for x in self._beam_parameterisations]
            params = [x for l in src_plists for x in l]
            global_p_list.extend(params)

        if self._xl_orientation_parameterisations:
            xlo_plists = [x.get_p() for x
                          in self._xl_orientation_parameterisations]
            params = [x for l in xlo_plists for x in l]
            global_p_list.extend(params)

        if self._xl_unit_cell_parameterisations:
            xluc_plists = [x.get_p() for x
                           in self._xl_unit_cell_parameterisations]
            params = [x for l in xluc_plists for x in l]
            global_p_list.extend(params)

        return global_p_list

    def set_p(self, vals):
        '''Set the parameter values of the contained models to the values in
        vals. This list must be of the same length as the result of get_p and
        must contain the parameter values in the same order! This order is to
        be maintained by any sensible refinement engine.'''

        assert len(vals) == len(self)
        it = iter(vals)

        if self._detector_parameterisations:
            for model in self._detector_parameterisations:
                tmp = [it.next() for i in range(model.num_free())]
                model.set_p(tmp)

        if self._beam_parameterisations:
            for model in self._beam_parameterisations:
                tmp = [it.next() for i in range(model.num_free())]
                model.set_p(tmp)

        if self._xl_orientation_parameterisations:
            for model in self._xl_orientation_parameterisations:
                tmp = [it.next() for i in range(model.num_free())]
                model.set_p(tmp)

        if self._xl_unit_cell_parameterisations:
            for model in self._xl_unit_cell_parameterisations:
                tmp = [it.next() for i in range(model.num_free())]
                model.set_p(tmp)

    def _prepare(self):
        '''Cache required quantities that are not dependent on hkl'''

        # Note, the reflection prediction code should be improved to also
        # provide detector + sensor numbers for each prediction, so that we can
        # have multiple sensor parameterisations and only calculate derivatives
        # for the one sensor that the ray does in fact intersect. We then need
        # a way to look up, from a detector + sensor number, which detector
        # parameterisation object refers to that sensor. Ideally this would be
        # done without requiring a search through all of
        # self._detector_parameterisations each time.
        #
        # For now, assume there is only one sensor, and it is parameterised by
        # the first entry in self._detector_parameterisations.

        ### Obtain various quantities of interest from the experimental model

        # Here we irrevocably choose the Panel that this reflection intersects,
        # currently hard-coding it to the first (only) Panel.
        self._D = matrix.sqr(self._detector[0].get_D_matrix())
        self._s0 = matrix.col(self._beam.get_s0())
        self._U = self._crystal.get_U()
        self._B = self._crystal.get_B()
        self._axis = matrix.col(self._gonio.get_rotation_axis())

    def get_gradients(self, h, s, phi):
        '''
        Calculate gradients of the prediction formula with respect to each
        of the parameters of the contained models, for the reflection with
        scattering vector s.

        To be implemented by a derived class, which determines the space of the
        prediction formula (e.g. we calculate dX/dp, dY/dp, dphi/dp for the
        prediction formula expressed in detector space, but components of
        d\vec{r}/dp for the prediction formula in reciprocal space
        '''

        self._prepare()

        return self._get_gradients_core(h, s, phi)

    def get_multi_gradients(self, match_list):
        '''
        perform the functionality of get_gradients but for a list of many
        reflections in one call in the form of a list of obs_pred_match objects
        (see target.py). The advantage of this is that _prepare needs only be
        called once.
        '''

        # This is effectively calculating the Jacobian and perhaps should be
        # renamed as such (and returned as a matrix not a list of lists)

        self._prepare()

        return [self._get_gradients_core(m.H, m.Sc, m.Phic) for m in match_list]

class detector_space_prediction_parameterisation_py(prediction_parameterisation):
    '''
    Concrete class that inherits functionaility of the
    prediction_parameterisation parent class and provides a detector space
    implementation of the get_gradients function.

    Not yet safe for multiple sensor detectors.

    Python version (slow)
    '''

    def _get_gradients_core(self, h, s, phi):

        '''Calculate gradients of the prediction formula with respect to each
        of the parameters of the contained models, for reflection h with
        scattering vector s that reflects at rotation angle phi. That is,
        calculate dX/dp, dY/dp and dphi/dp'''

        ### Calculate various quantities of interest for this reflection

        R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)

        # pv is the 'projection vector' for the reflection s.
        s = matrix.col(s)
        pv = self._D * s
        # r is the reciprocal lattice vector, in the lab frame
        r = R * self._U * self._B * h

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
            print "U =",self._U
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
            for idet, det in enumerate(self._detector_parameterisations):
                if idet == 0:
                    dd_ddet_p = det.get_ds_dp()
                    dpv_ddet_p = [- self._D * (dd_ddet_p[i]) * pv for i
                                     in range(len(dd_ddet_p))]
                else:
                    dpv_ddet_p = [matrix.col((0., 0., 0.))] * len(dd_ddet_p)

                dphi_ddet_p = [0.] * len(dd_ddet_p)

                dpv_dp.extend(dpv_ddet_p)
                dphi_dp.extend(dphi_ddet_p)

        # Calc derivatives of pv and phi wrt each parameter of each beam
        # parameterisation that is present.
        if self._beam_parameterisations:
            for src in self._beam_parameterisations:
                ds0_dsrc_p = src.get_ds_dp()
                dphi_dsrc_p = [- r.dot(ds0_dsrc_p[i]) / e_r_s0 for i
                                  in range(len(ds0_dsrc_p))]
                dpv_dsrc_p = [self._D * (e_X_r * dphi_dsrc_p[i] + ds0_dsrc_p[i]) for i in range(len(ds0_dsrc_p))]

                dpv_dp.extend(dpv_dsrc_p)
                dphi_dp.extend(dphi_dsrc_p)

        # Calc derivatives of pv and phi wrt each parameter of each crystal
        # orientation parameterisation that is present.
        if self._xl_orientation_parameterisations:
            for xlo in self._xl_orientation_parameterisations:
                dU_dxlo_p = xlo.get_ds_dp()

                dr_dxlo_p = [R * dU_dxlo_p[i] * self._B * h for i in range(len(dU_dxlo_p))]

                dphi_dxlo_p = [- der.dot(s) / e_r_s0 for der in dr_dxlo_p]

                dpv_dxlo_p = [self._D * (dr_dxlo_p[i] + e_X_r * dphi_dxlo_p[i]) for i in range(len(dphi_dxlo_p))]

                dpv_dp.extend(dpv_dxlo_p)
                dphi_dp.extend(dphi_dxlo_p)

        # Now derivatives of pv and phi wrt each parameter of each crystal unit
        # cell parameterisation that is present.
        if self._xl_unit_cell_parameterisations:
            for xluc in self._xl_unit_cell_parameterisations:
                dB_dxluc_p = xluc.get_ds_dp()

                dr_dxluc_p = [R * self._U * dB_dxluc_p[i] * h for i
                                  in range(len(dB_dxluc_p))]

                dphi_dxluc_p = [- der.dot(s) / e_r_s0 for der in dr_dxluc_p]

                dpv_dxluc_p = [self._D * (dr_dxluc_p[i] + e_X_r * dphi_dxluc_p[i]) for i in range(len(dr_dxluc_p))]

                dpv_dp.extend(dpv_dxluc_p)
                dphi_dp.extend(dphi_dxluc_p)

        # calculate positional derivatives from d[pv]/dp
        pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e) for e in dpv_dp]
        dX_dp, dY_dp = zip(*pos_grad)

        return zip(dX_dp, dY_dp, dphi_dp)

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

class detector_space_prediction_parameterisation(prediction_parameterisation):
    '''
    Concrete class that inherits functionaility of the
    prediction_parameterisation parent class and provides a detector space
    implementation of the get_gradients function.

    Not yet safe for multiple sensor detectors.
    '''

    def _get_gradients_core(self, h, s, phi):

        '''Calculate gradients of the prediction formula with respect to each
        of the parameters of the contained models, for reflection h with
        scattering vector s that reflects at rotation angle phi. That is,
        calculate dX/dp, dY/dp and dphi/dp'''

        from libtbx.test_utils import approx_equal

        ### Calculate various quantities of interest for this reflection

        R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)

        # pv is the 'projection vector' for the reflection s.
        s = matrix.col(s)
        pv = self._D * s
        # r is the reciprocal lattice vector, in the lab frame
        r = R * self._U * self._B * h

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
            print "U =",self._U
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
            for idet, det in enumerate(self._detector_parameterisations):
                if idet == 0:
                    dd_ddet_p = det.get_ds_dp()

                    #dpv_ddet_p = [- self._D * (dd_ddet_p[i]) * pv for i
                    #                 in range(len(dd_ddet_p))]

                    dd_ddet_p_temp = flex_mat3_double(len(dd_ddet_p))
                    for i, dd in enumerate(dd_ddet_p):
                        dd_ddet_p_temp[i] = dd.elems
                    dd_ddet_p = dd_ddet_p_temp

                    #print "detector_pv_derivative"
                    dpv_ddet_p = detector_pv_derivative(
                                    self._D, dd_ddet_p, pv)
                    #for d1, d2 in zip(dpv_ddet_p, dpv_ddet_p_2):
                    #    d11 = matrix.col(d1)
                    #    d22 = matrix.col(d2)
                    #    assert approx_equal(d11, d22, eps = 1.e-10)

                else:
                    dpv_ddet_p = [matrix.col((0., 0., 0.))] * len(dd_ddet_p)

                dphi_ddet_p = [0.] * len(dd_ddet_p)

                dpv_dp.extend(dpv_ddet_p)
                dphi_dp.extend(dphi_ddet_p)

        # Calc derivatives of pv and phi wrt each parameter of each beam
        # parameterisation that is present.
        if self._beam_parameterisations:
            for src in self._beam_parameterisations:
                ds0_dsrc_p = src.get_ds_dp()
                #dphi_dsrc_p = [- r.dot(ds0_dsrc_p[i]) / e_r_s0 for i
                #                  in range(len(ds0_dsrc_p))]
                #dpv_dsrc_p = [self._D * (e_X_r * dphi_dsrc_p[i] + ds0_dsrc_p[i]) for i in range(len(ds0_dsrc_p))]

                #print "beam_phi_derivative"
                dphi_dsrc_p = beam_phi_derivative(
                                    r, flex.vec3_double(ds0_dsrc_p), e_r_s0)

                #print "beam_pv_derivative"
                dpv_dsrc_p = beam_pv_derivative(
                                    self._D, e_X_r, dphi_dsrc_p,
                                    flex.vec3_double(ds0_dsrc_p))

                #for d1, d2 in zip(dphi_dsrc_p, dphi_dsrc_p_2):
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                #for d1, d2 in zip(dpv_dsrc_p, dpv_dsrc_p_2):
                #    d11 = matrix.col(d1)
                #    d22 = matrix.col(d2)
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                dpv_dp.extend(dpv_dsrc_p)
                dphi_dp.extend(dphi_dsrc_p)

        # Calc derivatives of pv and phi wrt each parameter of each crystal
        # orientation parameterisation that is present.
        if self._xl_orientation_parameterisations:
            for xlo in self._xl_orientation_parameterisations:
                dU_dxlo_p = xlo.get_ds_dp()

                #dr_dxlo_p = [R * dU_dxlo_p[i] * self._B * h for i in range(len(dU_dxlo_p))]

                #print "crystal_orientation_r_derivative"
                h2 = (int(h[0]), int(h[1]), int(h[2]))
                dr_dxlo_p = crystal_orientation_r_derivative(
                                R.elems, flex_mat3_double(dU_dxlo_p),
                                self._B.elems, h2)

                #dphi_dxlo_p = [- der.dot(s) / e_r_s0 for der in dr_dxlo_p]

                #print "crystal_orientation_phi_derivative"
                dphi_dxlo_p = crystal_orientation_phi_derivative(
                                    flex.vec3_double(dr_dxlo_p),
                                    s, e_r_s0)

                #dpv_dxlo_p = [self._D * (dr_dxlo_p[i] + e_X_r * dphi_dxlo_p[i]) for i in range(len(dphi_dxlo_p))]

                #print "crystal_orientation_pv_derivative"
                dpv_dxlo_p = crystal_orientation_pv_derivative(
                                    self._D.elems, dr_dxlo_p,
                                    e_X_r.elems, dphi_dxlo_p)

                #for d1, d2 in zip(dr_dxlo_p, dr_dxlo_p_2):
                #    d11 = matrix.col(d1)
                #    d22 = matrix.col(d2)
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                #for d1, d2 in zip(dphi_dxlo_p, dphi_dxlo_p_2):
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                #for d1, d2 in zip(dpv_dxlo_p, dpv_dxlo_p_2):
                #    d11 = matrix.col(d1)
                #    d22 = matrix.col(d2)
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                dpv_dp.extend(dpv_dxlo_p)
                dphi_dp.extend(dphi_dxlo_p)

        # Now derivatives of pv and phi wrt each parameter of each crystal unit
        # cell parameterisation that is present.
        if self._xl_unit_cell_parameterisations:
            for xluc in self._xl_unit_cell_parameterisations:
                dB_dxluc_p = xluc.get_ds_dp()

                #dr_dxluc_p = [R * self._U * dB_dxluc_p[i] * h for i
                #                  in range(len(dB_dxluc_p))]

                #print "crystal_cell_r_derivative"
                h2 = (int(h[0]), int(h[1]), int(h[2]))
                dr_dxluc_p = crystal_cell_r_derivative(
                                R.elems, self._U.elems,
                                flex_mat3_double(dB_dxluc_p), h2)

                #dphi_dxluc_p = [- der.dot(s) / e_r_s0 for der in dr_dxluc_p]

                #print "crystal_cell_phi_derivative"
                dphi_dxluc_p = crystal_cell_phi_derivative(
                                    dr_dxluc_p, s, e_r_s0)

                #dpv_dxluc_p = [self._D * (dr_dxluc_p[i] + e_X_r * dphi_dxluc_p[i]) for i in range(len(dr_dxluc_p))]

                #print "crystal_cell_pv_derivative"
                dpv_dxluc_p = crystal_cell_pv_derivative(
                                    self._D.elems, dr_dxluc_p,
                                    e_X_r.elems, dphi_dxluc_p)

                #for d1, d2 in zip(dr_dxluc_p, dr_dxluc_p_2):
                #    d11 = matrix.col(d1)
                #    d22 = matrix.col(d2)
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                #for d1, d2 in zip(dphi_dxluc_p, dphi_dxluc_p_2):
                #    assert approx_equal(d11, d22, eps = 1.e-10)

                #for d1, d2 in zip(dpv_dxluc_p, dpv_dxluc_p_2):
                #    d11 = matrix.col(d1)
                #    d22 = matrix.col(d2)
                #    assert approx_equal(d11, d22, eps = 1.e-10)


                dpv_dp.extend(dpv_dxluc_p)
                dphi_dp.extend(dphi_dxluc_p)

        # calculate positional derivatives from d[pv]/dp
        pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e) for e in dpv_dp]
        dX_dp, dY_dp = zip(*pos_grad)

        return zip(dX_dp, dY_dp, dphi_dp)

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
