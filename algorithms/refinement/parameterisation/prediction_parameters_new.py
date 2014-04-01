#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

#### Python and general cctbx imports

from __future__ import division
from scitbx import matrix

#### DIALS imports

from dials.array_family import flex
from dials_refinement_helpers_ext import *

"""This version of PredictionParameterisation vectorises calculations over
reflections, using flex arrays"""

class PredictionParameterisation(object):
  """
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

  * Detector parameterisation
  * Beam parameterisation
  * Crystal orientation parameterisation
  * Crystal unit cell parameterisation

  One of each must be supplied for each Experiment. These could be satisfied by
  a dummy class if no parameterisation is desired for some model.

  We also need access to the underlying models that are parameterised. The
  model parameterisation objects do not provide access to these models as it is
  not their job to do so. Instead we keep a separate reference to an
  ExperimentList that allows access to the relevant models.

  The goniometer is not yet parameterised, but we need it for the equations if
  we are doing parameterisation in X, Y, Phi space. Conversely, if
  parameterisation is only in X, Y space, the goniometer model is optional.

  A class implementing PredictionParameterisation is used by a Refinery
  object directly, which takes the list of parameters, and indirectly via a
  Target function object, which takes the list of derivatives and composes the
  derivatives of a Target function from them.

  """

  def __init__(self,
               experiments,
               detector_parameterisations = None,
               beam_parameterisations = None,
               xl_orientation_parameterisations = None,
               xl_unit_cell_parameterisations = None):

    # References to the underlying models
    self._experiments = experiments

    # Keep references to all parameterised models
    self._detector_parameterisations = detector_parameterisations
    self._beam_parameterisations = beam_parameterisations
    self._xl_orientation_parameterisations = \
        xl_orientation_parameterisations
    self._xl_unit_cell_parameterisations = \
        xl_unit_cell_parameterisations

    self._length = self._len()

    # Calculate Experiment to parameterisation mapping
    e2bp = dict([(ids, i) for i, dp in enumerate(beam_parameterisations) \
                 for ids in dp.get_experiment_ids()])
    e2xop = dict([(ids, i) for i, dp in enumerate(xl_orientation_parameterisations) \
                 for ids in dp.get_experiment_ids()])
    e2xucp = dict([(ids, i) for i, dp in enumerate(xl_unit_cell_parameterisations) \
                  for ids in dp.get_experiment_ids()])
    e2dp = dict([(ids, i) for i, dp in enumerate(detector_parameterisations) \
                 for ids in dp.get_experiment_ids()])
    from collections import namedtuple
    ParamSet = namedtuple('ParamSet', ['beam_param', 'xl_ori_param',
                                         'xl_uc_param', 'det_param'])
    self._exp_to_param = {i: ParamSet(e2bp[i], e2xop[i], e2xucp[i], e2dp[i]) \
                          for i, _ in enumerate(experiments)}

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

  def get_param_vals(self):
    """return a concatenated list of parameters from each of the components
    in the global model"""

    global_p_list = []
    if self._detector_parameterisations:
      det_plists = [x.get_param_vals() for x
                    in self._detector_parameterisations]
      params = [x for l in det_plists for x in l]
      global_p_list.extend(params)

    if self._beam_parameterisations:
      src_plists = [x.get_param_vals() for x
                    in self._beam_parameterisations]
      params = [x for l in src_plists for x in l]
      global_p_list.extend(params)

    if self._xl_orientation_parameterisations:
      xlo_plists = [x.get_param_vals() for x
                    in self._xl_orientation_parameterisations]
      params = [x for l in xlo_plists for x in l]
      global_p_list.extend(params)

    if self._xl_unit_cell_parameterisations:
      xluc_plists = [x.get_param_vals() for x
                     in self._xl_unit_cell_parameterisations]
      params = [x for l in xluc_plists for x in l]
      global_p_list.extend(params)

    return global_p_list

  def get_param_names(self):
    """Return a list of the names of parameters in the order they are
    concatenated. Useful for output to log files and debugging."""
    param_names = []
    if self._detector_parameterisations:
      det_param_name_lists = [x.get_param_names() for x in \
                         self._detector_parameterisations]
      names = ["Detector%d" % i + x for i, l \
               in enumerate(det_param_name_lists) for x in l]
      param_names.extend(names)

    if self._beam_parameterisations:
      beam_param_name_lists = [x.get_param_names() for x in \
                         self._beam_parameterisations]
      params = ["Beam%d" % i + x for i, l \
                in enumerate(beam_param_name_lists) for x in l]
      param_names.extend(params)

    if self._xl_orientation_parameterisations:
      xlo_param_name_lists = [x.get_param_names() for x
                    in self._xl_orientation_parameterisations]
      params = ["Crystal%d" % i + x for i, l \
                in enumerate(xlo_param_name_lists) for x in l]
      param_names.extend(params)

    if self._xl_unit_cell_parameterisations:
      xluc_param_name_lists = [x.get_param_names() for x
                     in self._xl_unit_cell_parameterisations]
      params = ["Crystal%d" % i + x for i, l \
                in enumerate(xluc_param_name_lists) for x in l]
      param_names.extend(params)

    return param_names

  def set_param_vals(self, vals):
    """Set the parameter values of the contained models to the values in
    vals. This list must be of the same length as the result of
    get_param_vals and must contain the parameter values in the same order!
    This order is to be maintained by any sensible refinement engine."""

    assert len(vals) == len(self)
    it = iter(vals)

    if self._detector_parameterisations:
      for model in self._detector_parameterisations:
        tmp = [it.next() for i in range(model.num_free())]
        model.set_param_vals(tmp)

    if self._beam_parameterisations:
      for model in self._beam_parameterisations:
        tmp = [it.next() for i in range(model.num_free())]
        model.set_param_vals(tmp)

    if self._xl_orientation_parameterisations:
      for model in self._xl_orientation_parameterisations:
        tmp = [it.next() for i in range(model.num_free())]
        model.set_param_vals(tmp)

    if self._xl_unit_cell_parameterisations:
      for model in self._xl_unit_cell_parameterisations:
        tmp = [it.next() for i in range(model.num_free())]
        model.set_param_vals(tmp)

  def get_gradients(self, reflections):
    """
    Calculate gradients of the prediction formula with respect to each
    of the parameters of the contained models, for all of the reflections.

    To be implemented by a derived class, which determines the space of the
    prediction formula (e.g. we calculate dX/dp, dY/dp, dphi/dp for the
    prediction formula expressed in detector space, but components of
    d\vec{r}/dp for the prediction formula in reciprocal space

    """

    ### Calculate various quantities of interest for the reflections

    # Set up arrays of values for each reflection
    # FIXME - completely replaces the prepare method...
    n = len(reflections)
    D = flex.mat3_double(n)
    s0 = flex.vec3_double(n)
    U = flex.mat3_double(n)
    B = flex.mat3_double(n)
    axis = flex.vec3_double(n)

    for iexp, exp in enumerate(self._experiments):

      sel = reflections['id'] == iexp
      isel = sel.iselection()

      # D matrix array
      panels = reflections['panel'].select(isel)
      for ipanel, D_mat in enumerate([p.get_D_matrix() for p in exp.detector]):
        subsel = isel.select(panels == ipanel)
        D.set_selected(subsel, D_mat)

      # s0 array
      s0.set_selected(isel, exp.beam.get_s0())

      # U array (scan varying version overrides this)
      U.set_selected(isel, exp.crystal.get_U())

      # B array (scan varying version overrides this)
      B.set_selected(isel, exp.crystal.get_B())

      # axis array
      if exp.goniometer:
        axis.set_selected(isel, exp.goniometer.get_rotation_axis())

    return self._get_gradients_core(reflections, D, s0, U, B, axis)

#  def get_gradients(self, h, s, phi, panel_id, obs_image_number=None,
#                    experiment_id=0):
#    """
#    Calculate gradients of the prediction formula with respect to each
#    of the parameters of the contained models, for the reflection with
#    scattering vector s that intersects panel with panel_id.
#
#    To be implemented by a derived class, which determines the space of the
#    prediction formula (e.g. we calculate dX/dp, dY/dp, dphi/dp for the
#    prediction formula expressed in detector space, but components of
#    d\vec{r}/dp for the prediction formula in reciprocal space
#
#    obs_image_number included to match the interface of a scan-
#    varying version of the class
#    """
#
#    # extract the right models
#    self._D = self._cache[experiment_id].D_mats[panel_id]
#    self._s0 = self._cache[experiment_id].s0
#    self._U = self._cache[experiment_id].U
#    self._B = self._cache[experiment_id].B
#    self._UB = self._cache[experiment_id].UB
#    self._axis = self._cache[experiment_id].axis
#
#    return self._get_gradients_core(h, s, phi, panel_id, experiment_id)


class XYPhiPredictionParameterisation(PredictionParameterisation):
  """
  Concrete class that inherits functionality of the
  PredictionParameterisation parent class and provides a detector space
  implementation of the get_gradients function.

  Untested for multiple sensor detectors.
  """

  def _get_gradients_core(self, reflections, D, s0, U, B, axis):
    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    that reflects at rotation angle phi with scattering vector s that
    intersects panel panel_id. That is, calculate dX/dp, dY/dp and
    dphi/dp"""

    # Spindle rotation matrices for every reflection
    #R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)
    #R = flex.mat3_double(len(reflections))
    # NB for now use flex.vec3_double.rotate_around_origin each time I need the
    # rotation matrix R.

    # pv is the 'projection vector' for the ray along s1.
    s1 = reflections['s1']
    pv = D * s1

    UB = U * B
    # r is the reciprocal lattice vector, in the lab frame
    h = reflections['miller_index'].as_vec3_double()
    phi_calc = reflections['xyzcal.mm'].parts()[2]
    r = (UB * h).rotate_around_origin(axis, phi_calc)

    # All of the derivatives of phi have a common denominator, given by
    # (e X r).s0, where e is the rotation axis. Calculate this once, here.
    e_X_r = axis.cross(r)
    e_r_s0 = (e_X_r).dot(s0)

    # Note that e_r_s0 -> 0 when the rotation axis, beam vector and
    # relp are coplanar. This occurs when a reflection just touches
    # the Ewald sphere.
    #
    # There is a relationship between e_r_s0 and zeta_factor.
    # Uncommenting the code below shows that
    # s0.(e X r) = zeta * |s X s0|

    #from dials.algorithms.reflection_basis import zeta_factor
    #from libtbx.test_utils import approx_equal
    #s = matrix.col(reflections['s1'][0])
    #z = zeta_factor(axis[0], s0[0], s)
    #ss0 = (s.cross(matrix.col(s0[0]))).length()
    #assert approx_equal(e_r_s0[0], z * ss0)

    # catch small values of e_r_s0
    e_r_s0_mag = flex.abs(e_r_s0)
    try:
      assert flex.min(e_r_s0_mag) > 1.e-6
    except AssertionError as e:
      imin = flex.min_index(e_r_s0_mag)
      print "(e X r).s0 too small:"
      print "for", (e_r_s0_mag <= 1.e-6).count(True), "reflections"
      print "out of", len(e_r_s0_mag), "total"
      print "such as", reflections['miller_index'][imin]
      print "with scattering vector", reflections['s1'][imin]
      print "where r =", r[imin]
      print "e =", axis[imin]
      print "s0 =", s0[imin]
      print ("this reflection forms angle with the equatorial plane "
             "normal:")
      vecn = matrix.col(s0[imin]).cross(matrix.col(axis[imin])).normalize()
      print matrix.col(reflections['s1'][imin]).accute_angle(vecn)
      raise e

    # Set up the lists of derivatives: a separate array over reflections for
    # each free parameter
    m = len(reflections)
    n = len(self) # number of free parameters
    dpv_dp = [flex.vec3_double(m, (0., 0., 0.)) for p in range(n)]
    dphi_dp = [flex.double(m, 0.) for p in range(n)]

    # set up return matrices FIXME might be better as 2D flex.grid, as commented out below
    # but for now use separate array for each free parameter
    #dX_dp = flex.double(flex.grid(m, n))
    #dY_dp = flex.double(flex.grid(m, n))
    #dphi_dp = flex.double(flex.grid(m, n))

    # loop over experiments
    for iexp, exp in enumerate(self._experiments):

      sel = reflections['id'] == iexp
      isel = sel.iselection()

      # identify which parameterisations to use for this experiment
      param_set = self._exp_to_param[iexp]
      beam_param_id = param_set.beam_param
      xl_ori_param_id = param_set.xl_ori_param
      xl_uc_param_id = param_set.xl_uc_param
      det_param_id = param_set.det_param

      # reset a pointer to the parameter number
      self._iparam = 0
      #print "for experiment", iexp, "the reflection indices are"
      #print list(isel)
      #print "the detector parameterisation id is", det_param_id
    #### FIXME from this point on

    ### Work through the parameterisations, calculating their contributions
    ### to derivatives d[pv]/dp and d[phi]/dp

      # Calculate derivatives of pv wrt each parameter of the detector
      # parameterisations. All derivatives of phi are zero for detector
      # parameters
      if self._detector_parameterisations:
        #self._detector_derivatives(dpv_dp, dphi_dp, pv, panel_id, det_param_id)
        self._detector_derivatives(reflections, isel, dpv_dp, D, pv, det_param_id, exp.detector)

      # Calc derivatives of pv and phi wrt each parameter of each beam
      # parameterisation that is present.
      if self._beam_parameterisations:
        #self._beam_derivatives(dpv_dp, dphi_dp, r, e_X_r, e_r_s0, beam_param_id)
        self._beam_derivatives(reflections, isel, dpv_dp, dphi_dp, r, e_X_r, e_r_s0, D, beam_param_id)

      # Calc derivatives of pv and phi wrt each parameter of each crystal
      # orientation parameterisation that is present.
      if self._xl_orientation_parameterisations:
        #self._xl_orientation_derivatives(dpv_dp, dphi_dp, R, h, s, \
        #                                   e_X_r, e_r_s0, xl_ori_param_id)
        self._xl_orientation_derivatives(reflections, isel, dpv_dp, dphi_dp, axis, phi_calc, h, s1, \
                                         e_X_r, e_r_s0, B, D, xl_ori_param_id, iexp)

      # Now derivatives of pv and phi wrt each parameter of each crystal unit
      # cell parameterisation that is present.
      if self._xl_unit_cell_parameterisations:
      #  self._xl_unit_cell_derivatives(dpv_dp, dphi_dp, R, h, s, \
      #                                   e_X_r, e_r_s0, xl_uc_param_id)
        self._xl_unit_cell_derivatives(reflections, isel, dpv_dp, dphi_dp, axis, phi_calc, h, s1, \
                                         e_X_r, e_r_s0, U, D, xl_uc_param_id, iexp)

      # calculate positional derivatives from d[pv]/dp
      dX_dp, dY_dp = self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, dpv_dp)
      #pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e)
      #            for e in dpv_dp]
      #dX_dp, dY_dp = zip(*pos_grad)

    return (dX_dp, dY_dp, dphi_dp)

  def _detector_derivatives(self, reflections, isel, dpv_dp, D, pv, det_param_id, detector):
    """helper function to extend the derivatives lists by
    derivatives of the detector parameterisations"""

    panels_this_exp = reflections['panel'].select(isel)

    # loop over all the detector parameterisations, even though we are only
    # setting values for one of them. We still need to move the _iparam pointer
    # for the others.
    for idp, dp in enumerate(self._detector_parameterisations):

      # Calculate gradients only for the correct detector parameterisation
      if idp == det_param_id:

        # loop through the panels in this detector
        for panel_id, panel in enumerate([p for p in detector]):

          # get the derivatives of detector d matrix for this panel
          dd_ddet_p = dp.get_ds_dp(multi_state_elt=panel_id)

          # get the right subset of array indices to set for this panel
          sub_isel = isel.select(panels_this_exp == panel_id)
          sub_pv = pv.select(sub_isel)
          sub_D = D.select(sub_isel)

          # loop through the parameters
          iparam = self._iparam
          for der in dd_ddet_p:

            # calculate the derivative of pv for this parameter
            dpv = (sub_D * (-1. * der).elems) * sub_pv

            # set values in the correct gradient array
            dpv_dp[iparam].set_selected(sub_isel, dpv)

            # increment the local parameter index pointer
            iparam += 1

        # increment the parameter index pointer to the last detector parameter
        self._iparam += dp.num_free()

      # For any other detector parameterisations, leave derivatives as zero
      else:

        # just increment the pointer
        self._iparam += dp.num_free()

    return

  def _beam_derivatives(self, reflections, isel, dpv_dp, dphi_dp, r, e_X_r, e_r_s0, D, beam_param_id):
    """helper function to extend the derivatives lists by
    derivatives of the beam parameterisations"""

    # loop over all the beam parameterisations, even though we are only setting
    # values for one of them. We still need to move the _iparam pointer for the
    # others.
    for ibp, bp in enumerate(self._beam_parameterisations):

      # Calculate gradients only for the correct beam parameterisation
      if ibp == beam_param_id:

        # get the derivatives of the beam vector wrt the parameters
        ds0_dbeam_p = bp.get_ds_dp()

        # select indices for the experiment of interest
        sub_r = r.select(isel)
        sub_e_X_r = e_X_r.select(isel)
        sub_e_r_s0 = e_r_s0.select(isel)
        sub_D = D.select(isel)

        # loop through the parameters
        for der in ds0_dbeam_p:

          # calculate the derivative of phi for this parameter
          dphi = (sub_r.dot(der.elems) / sub_e_r_s0) * -1.0

          # calculate the derivative of pv for this parameter
          dpv = sub_D * (sub_e_X_r * dphi + der)

          # set values in the correct gradient arrays
          dphi_dp[self._iparam].set_selected(isel, dphi)
          dpv_dp[self._iparam].set_selected(isel, dpv)

          # increment the parameter index pointer
          self._iparam += 1

      # For any other beam parameterisations, leave derivatives as zero
      else:

        # just increment the pointer
        self._iparam += bp.num_free()

    return

  def _xl_orientation_derivatives(self, reflections, isel, dpv_dp, dphi_dp, axis, phi_calc, h, s1, \
                                         e_X_r, e_r_s0, B, D, xl_ori_param_id, iexp):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations

    iexp is needed by the scan-varying override of this function only"""

    # loop over all the crystal orientation parameterisations, even though we
    # are only setting values for one of them. We still need to move the _iparam
    # pointer for the others.
    for ixlop, xlop in enumerate(self._xl_orientation_parameterisations):

      # Calculate gradients only for the correct xl orientation parameterisation
      if ixlop == xl_ori_param_id:

        # get derivatives of the U matrix wrt the parameters
        dU_dxlo_p = xlop.get_ds_dp()

        # select indices for the experiment of interest
        sub_h = h.select(isel)
        sub_s1 = s1.select(isel)
        #sub_r = r.select(isel)
        sub_e_X_r = e_X_r.select(isel)
        sub_e_r_s0 = e_r_s0.select(isel)
        sub_axis = axis.select(isel) # this is always the same vector for one
        # experiment! Seems wasteful to use an array
        sub_phi_calc = phi_calc.select(isel)
        sub_B = B.select(isel) # also always the same matrix for one experiment
        sub_D = D.select(isel)

        # loop through the parameters
        for der in dU_dxlo_p:

          der_mat = flex.mat3_double(len(sub_B), der.elems)
          # calculate the derivative of r for this parameter
          # FIXME COULD DO THIS BETTER WITH __rmul__?!
          tmp = der_mat * sub_B * sub_h
          dr = tmp.rotate_around_origin(sub_axis, sub_phi_calc)

          # calculate the derivative of phi for this parameter
          dphi = -1.0 * dr.dot(sub_s1) / sub_e_r_s0

          # calculate the derivative of pv for this parameter
          dpv = sub_D * (dr + sub_e_X_r * dphi)

          # set values in the correct gradient arrays
          dphi_dp[self._iparam].set_selected(isel, dphi)
          dpv_dp[self._iparam].set_selected(isel, dpv)

          # increment the parameter index pointer
          self._iparam += 1

      # For any other xl orientation parameterisations, leave derivatives as zero
      else:

        # just increment the pointer
        self._iparam += xlop.num_free()

    return

  def _xl_unit_cell_derivatives(self, reflections, isel, dpv_dp, dphi_dp, axis, phi_calc, h, s1, \
                                         e_X_r, e_r_s0, U, D, xl_uc_param_id, iexp):
    """helper function to extend the derivatives lists by
    derivatives of the crystal unit cell parameterisations

    iexp is needed by the scan-varying override of this function only"""

    for ixlucp, xlucp in enumerate(self._xl_unit_cell_parameterisations):

      # Calculate gradients only for the correct xl unit cell parameterisation
      if ixlucp == xl_uc_param_id:

        # get derivatives of the B matrix wrt the parameters
        dB_dxluc_p = xlucp.get_ds_dp()

        # select indices for the experiment of interest
        sub_h = h.select(isel)
        sub_s1 = s1.select(isel)
        #sub_r = r.select(isel)
        sub_e_X_r = e_X_r.select(isel)
        sub_e_r_s0 = e_r_s0.select(isel)
        sub_axis = axis.select(isel) # this is always the same vector for one
        # experiment! Seems wasteful to use an array
        sub_phi_calc = phi_calc.select(isel)
        sub_U = U.select(isel) # also always the same matrix for one experiment
        sub_D = D.select(isel)

        # loop through the parameters
        for der in dB_dxluc_p:

          der_mat = flex.mat3_double(len(sub_U), der.elems)
          # calculate the derivative of r for this parameter
          tmp = sub_U * der_mat * sub_h
          dr = tmp.rotate_around_origin(sub_axis, sub_phi_calc)

          # calculate the derivative of phi for this parameter
          dphi = -1.0 * dr.dot(sub_s1) / sub_e_r_s0

          # calculate the derivative of pv for this parameter
          dpv = sub_D * (dr + sub_e_X_r * dphi)

          # set values in the correct gradient arrays
          dphi_dp[self._iparam].set_selected(isel, dphi)
          dpv_dp[self._iparam].set_selected(isel, dpv)

          # increment the parameter index pointer
          self._iparam += 1

      # For any other xl unit cell parameterisations, leave derivatives as zero
      else:

        # just increment the pointer
        self._iparam += xlucp.num_free()

    return

  def _calc_dX_dp_and_dY_dp_from_dpv_dp(self, pv, dpv_dp):
    """helper function to calculate positional derivatives from
    dpv_dp using the quotient rule"""

    u, v, w = pv.parts()
    w2 = w**2

    dX_dp = []
    dY_dp = []

    for der in dpv_dp:
      du_dp, dv_dp, dw_dp = der.parts()

      dX_dp.append(du_dp / w - u * dw_dp / w2)
      dY_dp.append(dv_dp / w - v * dw_dp / w2)

    return dX_dp, dY_dp

