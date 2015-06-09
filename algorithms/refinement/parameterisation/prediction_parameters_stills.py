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

#### DIALS imports

from dials.array_family import flex
from dials.algorithms.refinement.parameterisation.prediction_parameters \
  import PredictionParameterisation, SparseGradientVectorMixin

#class StillsPredictionParameterisation(PredictionParameterisation):
#  """
#  Concrete class that inherits functionality of the PredictionParameterisation
#  parent class and provides a detector space implementation of the get_gradients
#  function for still images.
#
#  """
#
#  def _get_gradients_core(self, reflections, D, s0, U, B, axis):
#    """Calculate gradients of the prediction formula with respect to
#    each of the parameters of the contained models, for reflection h
#    with scattering vector s that intersects panel panel_id. That is,
#    calculate dX/dp, dY/dp and dDeltaPsi/dp. Ignore axis because these
#    are stills"""
#
#    # pv is the 'projection vector' for the ray along s1.
#    self._D = D
#    self._s1 = reflections['s1']
#    self._pv = D * self._s1
#
#    # q is the reciprocal lattice vector, in the lab frame
#    self._UB = U * B
#    self._U = U
#    self._B = B
#    self._h = reflections['miller_index'].as_vec3_double()
#    self._q = (self._UB * self._h)
#
#    # r is the reciprocal lattice vector rotated to the Ewald sphere
#    self._s0 = s0
#    self._r = self._s1 - self._s0
#
#    # we also need the unit directions q0 and s0u
#    self._q0 = self._q.each_normalize()
#    self._s0u = self._s0.each_normalize()
#
#    # e1 is the unit vector about which DeltaPsi rotation is defined
#    self._e1 = self._q0.cross(self._s0u).each_normalize()
#
#    # q1 completes an orthonormal set with q0 and e1
#    self._q1 = self._q0.cross(self._e1).each_normalize()
#
#    # c0 completes an orthonormal set with s0u and e1
#    self._c0 = self._s0u.cross(self._e1).each_normalize()
#
#    # calculate 'a' and 'b', the components of r in the directions -s0u and
#    # c0 respectively
#    self._a = -1.0 * self._r.dot(self._s0u)
#    self._b = self._r.dot(self._c0)
#
#    # we want the wavelength
#    self._wavelength = 1. / self._s0.norms()
#
#    # calculate chi and upsilon. In the paper these are called 'x' and 'y', but
#    # we rename here to avoid confusion with detector coordinates X and Y.
#    self._upsilon = self._r.dot(self._q1)
#    self._chi = self._r.dot(self._q0)
#
#    # Set up the lists of derivatives: a separate array over reflections for
#    # each free parameter
#    m = len(reflections)
#    n = len(self) # number of free parameters
#    dpv_dp = [flex.vec3_double(m, (0., 0., 0.)) for p in range(n)]
#    dchi_dp = [flex.double(m, 0.) for p in range(n)]
#    dupsilon_dp = [flex.double(m, 0.) for p in range(n)]
#
#    # loop over experiments
#    for iexp, exp in enumerate(self._experiments):
#
#      sel = reflections['id'] == iexp
#      isel = sel.iselection()
#
#      # identify which parameterisations to use for this experiment
#      param_set = self._exp_to_param[iexp]
#      beam_param_id = param_set.beam_param
#      xl_ori_param_id = param_set.xl_ori_param
#      xl_uc_param_id = param_set.xl_uc_param
#      det_param_id = param_set.det_param
#
#      # reset a pointer to the parameter number
#      self._iparam = 0
#
#    ### Work through the parameterisations, calculating their contributions
#    ### to derivatives d[pv]/dp and d[DeltaPsi]/dp
#
#      # Calculate derivatives of pv wrt each parameter of the detector
#      # parameterisations. All derivatives of DeltaPsi are zero for detector
#      # parameters
#      if self._detector_parameterisations:
#        self._detector_derivatives(reflections, isel, dpv_dp,
#                                   det_param_id, exp.detector)
#
#      # Calc derivatives of pv and DeltaPsi wrt each parameter of each beam
#      # parameterisation that is present.
#      if self._beam_parameterisations:
#        self._beam_derivatives(reflections, isel, dpv_dp, dchi_dp, dupsilon_dp,
#                               beam_param_id)
#
#      # Calc derivatives of pv and phi wrt each parameter of each crystal
#      # orientation parameterisation that is present.
#      if self._xl_orientation_parameterisations:
#        self._xl_orientation_derivatives(reflections, isel, dpv_dp,
#                                         dchi_dp, dupsilon_dp, xl_ori_param_id)
#
#      # Now derivatives of pv and phi wrt each parameter of each crystal unit
#      # cell parameterisation that is present.
#      if self._xl_unit_cell_parameterisations:
#        self._xl_unit_cell_derivatives(reflections, isel, dpv_dp,
#                                       dchi_dp, dupsilon_dp, xl_uc_param_id)
#
#      # calculate positional derivatives from d[pv]/dp
#      dX_dp, dY_dp = self._calc_dX_dp_and_dY_dp_from_dpv_dp(self._pv, dpv_dp)
#
#      # calculate angular derivatives from d[chi]/dp and d[upsilon]/dp
#      dDeltaPsi_dp = self._calc_dDeltaPsi_dp(self._chi, self._upsilon,
#                                             dchi_dp, dupsilon_dp)
#
#    return (dX_dp, dY_dp, dDeltaPsi_dp)
#
#
#  def _detector_derivatives(self, reflections, isel, dpv_dp,
#                            det_param_id, detector):
#    """helper function to extend the derivatives lists by derivatives of the
#    detector parameterisations"""
#
#    panels_this_exp = reflections['panel'].select(isel)
#
#    # loop over all the detector parameterisations, even though we are only
#    # setting values for one of them. We still need to move the _iparam pointer
#    # for the others.
#    for idp, dp in enumerate(self._detector_parameterisations):
#
#      # Calculate gradients only for the correct detector parameterisation
#      if idp == det_param_id:
#
#        # loop through the panels in this detector
#        for panel_id, panel in enumerate([p for p in detector]):
#
#          # get the derivatives of detector d matrix for this panel
#          dd_ddet_p = dp.get_ds_dp(multi_state_elt=panel_id)
#
#          # get the right subset of array indices to set for this panel
#          sub_isel = isel.select(panels_this_exp == panel_id)
#          pv = self._pv.select(sub_isel)
#          D = self._D.select(sub_isel)
#
#          # loop through the parameters
#          iparam = self._iparam
#          for der in dd_ddet_p:
#
#            # calculate the derivative of pv for this parameter
#            dpv = (D * (-1. * der).elems) * pv
#
#            # set values in the correct gradient array
#            dpv_dp[iparam].set_selected(sub_isel, dpv)
#
#            # increment the local parameter index pointer
#            iparam += 1
#
#        # increment the parameter index pointer to the last detector parameter
#        self._iparam += dp.num_free()
#
#      # For any other detector parameterisations, leave derivatives as zero
#      else:
#
#        # just increment the pointer
#        self._iparam += dp.num_free()
#
#    return
#
#  def _beam_derivatives(self, reflections, isel, dpv_dp,
#                        dchi_dp, dupsilon_dp, beam_param_id):
#    """helper function to extend the derivatives lists by derivatives of the
#    beam parameterisations"""
#
#    # loop over all the beam parameterisations, even though we are only setting
#    # values for one of them. We still need to move the _iparam pointer for the
#    # others.
#    for ibp, bp in enumerate(self._beam_parameterisations):
#
#      # Calculate gradients only for the correct beam parameterisation
#      if ibp == beam_param_id:
#
#        # get the derivatives of the beam vector wrt the parameters
#        ds0_dbeam_p = bp.get_ds_dp()
#
#        # select indices for the experiment of interest
#        s0u = self._s0u.select(isel)
#        wl = self._wavelength.select(isel)
#        q0 = self._q0.select(isel)
#        e1 = self._e1.select(isel)
#        a = self._a.select(isel)
#        b = self._b.select(isel)
#        q1 = self._q1.select(isel)
#        r = self._r.select(isel)
#        chi = self._chi.select(isel)
#        upsilon = self._upsilon.select(isel)
#        D = self._D.select(isel)
#
#        # loop through the parameters
#        for der in ds0_dbeam_p:
#
#          # repeat the derivative in an array
#          ds0 = flex.vec3_double(len(s0u), der.elems)
#
#          # we need the derivative of the unit beam direction too
#          ds0u = ds0 * wl
#
#          # calculate the derivatives of upsilon and chi for this parameter
#          de1 = q0.cross(ds0u)
#          dc0 = ds0u.cross(e1) + s0u.cross(de1)
#          dr = b * dc0 - a * ds0u
#          dq1 = q0.cross(de1)
#
#          dupsilon = dr.dot(q1) + r.dot(dq1)
#          dchi = dr.dot(q0)
#
#          # calculate the derivative of pv for this parameter
#          dpv = D * (ds0 + dr)
#
#          # set values in the correct gradient arrays
#          dchi_dp[self._iparam].set_selected(isel, dchi)
#          dupsilon_dp[self._iparam].set_selected(isel, dupsilon)
#          dpv_dp[self._iparam].set_selected(isel, dpv)
#
#          # increment the parameter index pointer
#          self._iparam += 1
#
#      # For any other beam parameterisations, leave derivatives as zero
#      else:
#
#        # just increment the pointer
#        self._iparam += bp.num_free()
#
#    return
#
#  def _xl_orientation_derivatives(self, reflections, isel, dpv_dp,
#                                  dchi_dp, dupsilon_dp, xl_ori_param_id):
#    """helper function to extend the derivatives lists by
#    derivatives of the crystal orientation parameterisations"""
#
#    # loop over all the crystal orientation parameterisations, even though we
#    # are only setting values for one of them. We still need to move the _iparam
#    # pointer for the others.
#    for ixlop, xlop in enumerate(self._xl_orientation_parameterisations):
#
#      # Calculate gradients only for the correct xl orientation parameterisation
#      if ixlop == xl_ori_param_id:
#
#        # get derivatives of the U matrix wrt the parameters
#        dU_dxlo_p = xlop.get_ds_dp()
#
#        # select indices for the experiment of interest
#        h = self._h.select(isel)
#        B = self._B.select(isel)
#        q = self._q.select(isel)
#        q0 = self._q0.select(isel)
#        s0u = self._s0u.select(isel)
#        wl = self._wavelength.select(isel)
#        b = self._b.select(isel)
#        c0 = self._c0.select(isel)
#        e1 = self._e1.select(isel)
#        q1 = self._q1.select(isel)
#        r = self._r.select(isel)
#        D = self._D.select(isel)
#
#        # loop through the parameters
#        for der in dU_dxlo_p:
#
#          der_mat = flex.mat3_double(len(B), der.elems)
#          # calculate the derivative of q for this parameter
#          dq = der_mat * B * h
#
#          # calculate the derivatives of upsilon and chi for this parameter
#          q_dot_dq = q.dot(dq)
#          dqq = 2.0 * q_dot_dq
#          q_scalar = q.norms()
#          qq = q_scalar * q_scalar
#          dq_scalar = dqq / q_scalar
#          dq0 = (q_scalar * dq - (q_dot_dq * q0)) / qq
#          de1 = dq0.cross(s0u)
#          dc0 = s0u.cross(de1)
#          da = wl * q_dot_dq
#          db = (2.0 - qq * wl * wl) * q_dot_dq / (2.0 * b)
#          dr = db * c0 + b * dc0 - da * s0u
#          dq1 = dq0.cross(e1) + q0.cross(de1)
#
#          dupsilon = dr.dot(q1) + r.dot(dq1)
#          dchi = dr.dot(q0) + r.dot(dq0)
#
#          # calculate the derivative of pv for this parameter
#          dpv = D * dq
#
#          # set values in the correct gradient arrays
#          dchi_dp[self._iparam].set_selected(isel, dchi)
#          dupsilon_dp[self._iparam].set_selected(isel, dupsilon)
#          dpv_dp[self._iparam].set_selected(isel, dpv)
#
#          # increment the parameter index pointer
#          self._iparam += 1
#
#      # For any other xl orientation parameterisations, leave derivatives as zero
#      else:
#
#        # just increment the pointer
#        self._iparam += xlop.num_free()
#
#    return
#
#  def _xl_unit_cell_derivatives(self, reflections, isel, dpv_dp,
#                                dchi_dp, dupsilon_dp, xl_uc_param_id):
#    """helper function to extend the derivatives lists by
#    derivatives of the crystal unit cell parameterisations"""
#
#
#    for ixlucp, xlucp in enumerate(self._xl_unit_cell_parameterisations):
#
#      # Calculate gradients only for the correct xl unit cell parameterisation
#      if ixlucp == xl_uc_param_id:
#
#        # get derivatives of the B matrix wrt the parameters
#        dB_dxluc_p = xlucp.get_ds_dp()
#
#        # select indices for the experiment of interest
#        h = self._h.select(isel)
#        U = self._U.select(isel)
#        q = self._q.select(isel)
#        q0 = self._q0.select(isel)
#        s0u = self._s0u.select(isel)
#        wl = self._wavelength.select(isel)
#        b = self._b.select(isel)
#        c0 = self._c0.select(isel)
#        e1 = self._e1.select(isel)
#        q1 = self._q1.select(isel)
#        r = self._r.select(isel)
#        D = self._D.select(isel)
#
#        # loop through the parameters
#        for der in dB_dxluc_p:
#
#          der_mat = flex.mat3_double(len(U), der.elems)
#          # calculate the derivative of q for this parameter
#          dq = U * der_mat * h
#
#          # calculate the derivatives of upsilon and chi for this parameter
#          q_dot_dq = q.dot(dq)
#          dqq = 2.0 * q_dot_dq
#          q_scalar = q.norms()
#          qq = q_scalar * q_scalar
#          dq_scalar = dqq / q_scalar
#          dq0 = (q_scalar * dq - (q_dot_dq * q0)) / qq
#          de1 = dq0.cross(s0u)
#          dc0 = s0u.cross(de1)
#          da = wl * q_dot_dq
#          db = (2.0 - qq * wl * wl) * q_dot_dq / (2.0 * b)
#          dr = db * c0 + b * dc0 - da * s0u
#          dq1 = dq0.cross(e1) + q0.cross(de1)
#
#          dupsilon = dr.dot(q1) + r.dot(dq1)
#          dchi = dr.dot(q0) + r.dot(dq0)
#
#          # calculate the derivative of pv for this parameter
#          dpv = D * dq
#
#          # set values in the correct gradient arrays
#          dchi_dp[self._iparam].set_selected(isel, dchi)
#          dupsilon_dp[self._iparam].set_selected(isel, dupsilon)
#          dpv_dp[self._iparam].set_selected(isel, dpv)
#
#          # increment the parameter index pointer
#          self._iparam += 1
#
#      # For any other xl unit cell parameterisations, leave derivatives as zero
#      else:
#
#        # just increment the pointer
#        self._iparam += xlucp.num_free()
#
#    return
#
#  @staticmethod
#  def _calc_dX_dp_and_dY_dp_from_dpv_dp(pv, dpv_dp):
#    """helper function to calculate positional derivatives from
#    dpv_dp using the quotient rule"""
#
#    u, v, w = pv.parts()
#
#    # precalculate for efficiency
#    w_inv = 1/w
#    u_w_inv = u * w_inv
#    v_w_inv = v * w_inv
#
#    dX_dp = []
#    dY_dp = []
#
#    for der in dpv_dp:
#      du_dp, dv_dp, dw_dp = der.parts()
#
#      dX_dp.append(w_inv * (du_dp - dw_dp * u_w_inv))
#      dY_dp.append(w_inv * (dv_dp - dw_dp * v_w_inv))
#
#    return dX_dp, dY_dp
#
#  @staticmethod
#  def _calc_dDeltaPsi_dp(chi, upsilon, dchi_dp, dupsilon_dp):
#    """helper function to calculate DeltaPsi derivatives from
#    d[chi]/dp and d[upsilon]/dp"""
#
#    # Using DeltaPsi = -arctan(upsilon/chi)
#
#    # precalculate shared coefficient
#    c = 1.0/(chi * chi + upsilon * upsilon)
#
#    dDeltaPsi_dp = []
#
#    for dchi, dupsilon in zip(dchi_dp, dupsilon_dp):
#
#      dDeltaPsi_dp.append(c * (upsilon * dchi - chi * dupsilon))
#
#    return dDeltaPsi_dp


class StillsPredictionParameterisation(PredictionParameterisation):
  """
  Concrete class that inherits functionality of the PredictionParameterisation
  parent class and provides a detector space implementation of the get_gradients
  function for still images.

  """

  _grad_names = ("dX_dp", "dY_dp", "dDeltaPsi_dp")

  def _get_gradients_core(self, reflections, D, s0, U, B, axis, fixed_rotation, callback=None):
    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    with scattering vector s that intersects panel panel_id. That is,
    calculate dX/dp, dY/dp and dDeltaPsi/dp. Ignore axis and fixed_rotation
    because these are stills"""

    # pv is the 'projection vector' for the ray along s1.
    self._D = D
    self._s1 = reflections['s1']
    self._pv = D * self._s1

    # also need quantities derived from pv, precalculated for efficiency
    u, v, w = self._pv.parts()
    self._w_inv = 1/w
    self._u_w_inv = u * self._w_inv
    self._v_w_inv = v * self._w_inv

    self._DeltaPsi = reflections['delpsical.rad']

    # q is the reciprocal lattice vector, in the lab frame
    self._UB = U * B
    self._U = U
    self._B = B
    self._h = reflections['miller_index'].as_vec3_double()
    self._q = (self._UB * self._h)
    self._q_scalar = self._q.norms()
    self._qq = self._q_scalar * self._q_scalar

    # r is the reciprocal lattice vector rotated to the Ewald sphere
    self._s0 = s0
    self._r = self._s1 - self._s0

    # we also need the unit directions q0 and s0u
    self._q0 = self._q.each_normalize()
    self._s0u = self._s0.each_normalize()

    # e1 is the unit vector about which DeltaPsi rotation is defined
    self._e1 = self._q0.cross(self._s0u).each_normalize()

    # q1 completes an orthonormal set with q0 and e1
    self._q1 = self._q0.cross(self._e1).each_normalize()

    # we want the wavelength
    self._wavelength = 1. / self._s0.norms()

    # Set up empty list in which to store gradients
    m = len(reflections)
    results = []

    # determine experiment to indices mappings once, here
    experiment_to_idx = []
    for iexp, exp in enumerate(self._experiments):

      sel = reflections['id'] == iexp
      isel = sel.iselection()
      experiment_to_idx.append(isel)

    # reset a pointer to the parameter number
    self._iparam = 0

    ### Work through the parameterisations, calculating their contributions
    ### to derivatives d[pv]/dp and d[DeltaPsi]/dp

    # loop over the detector parameterisations
    for dp in self._detector_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in dp.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Access the detector model being parameterised
      detector = dp.get_model()

      # Get panel numbers of the affected reflections
      panel = reflections['panel'].select(isel)

      # Extend derivative vectors for this detector parameterisation
      results = self._extend_gradient_vectors(results, m, dp.num_free(),
        keys=self._grad_names)

      # loop through the panels in this detector
      for panel_id, _ in enumerate(detector):

        # get the right subset of array indices to set for this panel
        sub_isel = isel.select(panel == panel_id)
        if len(sub_isel) == 0:
          # if no reflections intersect this panel, skip calculation
          continue
        sub_pv = self._pv.select(sub_isel)
        sub_D = self._D.select(sub_isel)
        dpv_ddet_p = self._detector_derivatives(dp, sub_pv, sub_D, panel_id)

        # convert to dX/dp, dY/dp and assign the elements of the vectors
        # corresponding to this experiment and panel
        sub_w_inv = self._w_inv.select(sub_isel)
        sub_u_w_inv = self._u_w_inv.select(sub_isel)
        sub_v_w_inv = self._v_w_inv.select(sub_isel)
        dX_ddet_p, dY_ddet_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
          sub_w_inv, sub_u_w_inv, sub_v_w_inv, dpv_ddet_p)

        # use a local parameter index pointer because we set all derivatives
        # for this panel before moving on to the next
        iparam = self._iparam
        for dX, dY in zip(dX_ddet_p, dY_ddet_p):
          results[iparam][self._grad_names[0]].set_selected(sub_isel, dX)
          results[iparam][self._grad_names[1]].set_selected(sub_isel, dY)
          # increment the local parameter index pointer
          iparam += 1

      if callback is not None:
        iparam = self._iparam
        for i in range(dp.num_free()):
          results[iparam] = callback(results[iparam])
          iparam += 1

      # increment the parameter index pointer to the last detector parameter
      self._iparam += dp.num_free()

    # loop over the beam parameterisations
    for bp in self._beam_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in bp.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Extend derivative vectors for this beam parameterisation
      results = self._extend_gradient_vectors(results, m, bp.num_free(),
        keys=self._grad_names)

      if len(isel) == 0:
        # if no reflections are in this experiment, skip calculation
        self._iparam += bp.num_free()
        continue

      # Get required data from those reflections
      s0 = self._s0.select(isel)
      s0u = self._s0u.select(isel)
      wl = self._wavelength.select(isel)
      r = self._r.select(isel)
      e1 = self._e1.select(isel)
      q = self._q.select(isel)
      q0 = self._q0.select(isel)
      DeltaPsi = self._DeltaPsi.select(isel)
      D = self._D.select(isel)

      w_inv = self._w_inv.select(isel)
      u_w_inv = self._u_w_inv.select(isel)
      v_w_inv = self._v_w_inv.select(isel)

      dpv_dbeam_p, ddelpsi_dbeam_p = self._beam_derivatives(
          bp, s0, s0u, wl, r, e1, q, q0, DeltaPsi, D)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this parameterisation
      dX_dbeam_p, dY_dbeam_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dbeam_p)
      for dX, dY, dDeltaPsi in zip(dX_dbeam_p, dY_dbeam_p, ddelpsi_dbeam_p):
        results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        results[self._iparam][self._grad_names[2]].set_selected(isel, dDeltaPsi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    # loop over the crystal orientation parameterisations
    for xlop in self._xl_orientation_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in xlop.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Extend derivative vectors for this crystal orientation parameterisation
      results = self._extend_gradient_vectors(results, m, xlop.num_free(),
        keys=self._grad_names)

      if len(isel) == 0:
        # if no reflections are in this experiment, skip calculation
        self._iparam += xlop.num_free()
        continue

      # Get required data from those reflections
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

      w_inv = self._w_inv.select(isel)
      u_w_inv = self._u_w_inv.select(isel)
      v_w_inv = self._v_w_inv.select(isel)

      dpv_dxlo_p, ddelpsi_dxlo_p = self._xl_orientation_derivatives(
        xlop, B, h, e1, DeltaPsi, s1, q, q_scalar, qq, q0, r, s0, s0u, D)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dxlo_p, dY_dxlo_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dxlo_p)
      for dX, dY, dDeltaPsi in zip(dX_dxlo_p, dY_dxlo_p, ddelpsi_dxlo_p):
        results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        results[self._iparam][self._grad_names[2]].set_selected(isel, dDeltaPsi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    # loop over the crystal unit cell parameterisations
    for xlucp in self._xl_unit_cell_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in xlucp.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Extend derivative vectors for this crystal unit cell parameterisation
      results = self._extend_gradient_vectors(results, m, xlucp.num_free(),
        keys=self._grad_names)

      if len(isel) == 0:
        # if no reflections are in this experiment, skip calculation
        self._iparam += xlucp.num_free()
        continue

      # Get required data from those reflections
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

      w_inv = self._w_inv.select(isel)
      u_w_inv = self._u_w_inv.select(isel)
      v_w_inv = self._v_w_inv.select(isel)

      dpv_dxluc_p, ddelpsi_dxluc_p =  self._xl_unit_cell_derivatives(
        xlucp, U, h, e1, DeltaPsi, s1, q, q_scalar, qq, q0, r, s0, s0u, D)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dxluc_p, dY_dxluc_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dxluc_p)
      for dX, dY, dDeltaPsi in zip(dX_dxluc_p, dY_dxluc_p, ddelpsi_dxluc_p):
        results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        results[self._iparam][self._grad_names[2]].set_selected(isel, dDeltaPsi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    return results


  def _detector_derivatives(self, dp, pv, D, panel_id):
    """helper function to convert derivatives of the detector state to
    derivatives of the vector pv"""

    # get the derivatives of detector d matrix for this panel
    dd_ddet_p = dp.get_ds_dp(multi_state_elt=panel_id)

    # calculate the derivative of pv for this parameter
    dpv_ddet_p = [(D * (-1. * der).elems) * pv for der in dd_ddet_p]

    return dpv_ddet_p

  def _beam_derivatives(self, bp, s0, s0u, wl, r, e1, q, q0, DeltaPsi, D):
    """helper function to extend the derivatives lists by derivatives of the
    beam parameterisations"""

    # get the derivatives of the beam vector wrt the parameters
    ds0_dbeam_p = bp.get_ds_dp()

    dDeltaPsi_dp = []
    dpv_dp = []

    # loop through the parameters
    for der in ds0_dbeam_p:

      # repeat the derivative in an array
      ds0 = flex.vec3_double(len(s0u), der.elems)

      # we need the derivative of the unit beam direction too
      ds0u = ds0 * wl

      # calculate the derivative of DeltaPsi for this parameter
      dDeltaPsi = -1.0 * (r.dot(ds0)) / (e1.cross(r).dot(s0))
      dDeltaPsi_dp.append(dDeltaPsi)

      # calculate the partial derivative of r wrt change in rotation
      # axis e1 by finite differences
      de1 = q0.cross(ds0u)
      dp = 1.e-8 # finite step size for the parameter
      del_e1 = de1 * dp
      e1f = e1 + del_e1 * 0.5
      rfwd = q.rotate_around_origin(e1f, DeltaPsi)
      e1r = e1 - del_e1 * 0.5
      rrev = q.rotate_around_origin(e1r, DeltaPsi)
      drde_dedp = (rfwd - rrev) * (1 / dp)

      # calculate the derivative of pv for this parameter
      dpv_dp.append(D * (ds0 + e1.cross(r) * dDeltaPsi + drde_dedp))

    return dpv_dp, dDeltaPsi_dp

  def _xl_orientation_derivatives(self, xlop, B, h, e1, DeltaPsi, s1, q,
                                  q_scalar, qq, q0, r, s0, s0u, D):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations"""

    # get derivatives of the U matrix wrt the parameters
    dU_dxlo_p = xlop.get_ds_dp()

    dDeltaPsi_dp = []
    dpv_dp = []

    # loop through the parameters
    for der in dU_dxlo_p:

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
      de1 = dq0.cross(s0u)

      # calculate the partial derivative of r wrt change in rotation
      # axis e1 by finite differences
      dp = 1.e-8
      del_e1 = de1 * dp
      e1f = e1 + del_e1 * 0.5
      rfwd = q.rotate_around_origin(e1f, DeltaPsi)
      e1r = e1 - del_e1 * 0.5
      rrev = q.rotate_around_origin(e1r, DeltaPsi)
      drde_dedp = (rfwd - rrev) * (1 / dp)

      # calculate the derivative of pv for this parameter
      dpv_dp.append(D * (dr + e1.cross(r) * dDeltaPsi + drde_dedp))

    return dpv_dp, dDeltaPsi_dp

  def _xl_unit_cell_derivatives(self, xlucp, U, h, e1, DeltaPsi, s1, q,
                                  q_scalar, qq, q0, r, s0, s0u, D):
    """helper function to extend the derivatives lists by
    derivatives of the crystal unit cell parameterisations"""

    # get derivatives of the B matrix wrt the parameters
    dB_dxluc_p = xlucp.get_ds_dp()

    dDeltaPsi_dp = []
    dpv_dp = []

    # loop through the parameters
    for der in dB_dxluc_p:

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
      de1 = dq0.cross(s0u)

      # calculate the partial derivative of r wrt change in rotation
      # axis e1 by finite differences
      dp = 1.e-8
      del_e1 = de1 * dp
      e1f = e1 + del_e1 * 0.5
      rfwd = q.rotate_around_origin(e1f , DeltaPsi)
      e1r = e1 - del_e1 * 0.5
      rrev = q.rotate_around_origin(e1r , DeltaPsi)
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
      du_dp, dv_dp, dw_dp = der.parts()

      dX_dp.append(w_inv * (du_dp - dw_dp * u_w_inv))
      dY_dp.append(w_inv * (dv_dp - dw_dp * v_w_inv))

    return dX_dp, dY_dp


class StillsPredictionParameterisationSparse(SparseGradientVectorMixin,
  StillsPredictionParameterisation):
  """A version of StillsPredictionParameterisation that uses a sparse matrix
  data structure for memory efficiency when there are a large number of
  Experiments"""

  pass
