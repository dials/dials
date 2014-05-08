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
from dials.algorithms.refinement.parameterisation.prediction_parameters \
  import PredictionParameterisation

class StillsPredictionParameterisation(PredictionParameterisation):
  """
  Concrete class that inherits functionality of the PredictionParameterisation
  parent class and provides a detector space implementation of the get_gradients
  function for still images.

  """

  def _get_gradients_core(self, reflections, D, s0, U, B, axis):
    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    with scattering vector s that intersects panel panel_id. That is,
    calculate dX/dp, dY/dp and dDeltaPsi/dp. Ignore axis because these
    are stills"""

    # pv is the 'projection vector' for the ray along s1.
    s1 = reflections['s1']
    pv = D * s1

    # q is the reciprocal lattice vector, in the lab frame
    h = reflections['miller_index'].as_vec3_double()
    UB = U * B
    q = (UB * h)

    # r is the reciprocal lattice vector rotated to the Ewald sphere
    r = s1 - s0

    # we also need the unit directions q0 and s0u
    q0 = q.each_normalize()
    s0u = s0.each_normalize()

    # e1 is the unit vector about which DeltaPsi rotation is defined
    e1 = q0.cross(s0u).each_normalize()

    # q1 completes an orthonormal set with q0 and e1
    q1 = q0.cross(e1).each_normalize()

    # c0 completes an orthonormal set with s0u and e1
    c0 = s0u.cross(e1).each_normalize()

    # calculate 'a' and 'b', the components of r in the directions -s0u and
    # c0 respectively
    a = -1.0 * r.dot(s0u)
    b = r.dot(c0)

    # Set up the lists of derivatives: a separate array over reflections for
    # each free parameter
    m = len(reflections)
    n = len(self) # number of free parameters
    dpv_dp = [flex.vec3_double(m, (0., 0., 0.)) for p in range(n)]
    dDeltaPsi_dp = [flex.double(m, 0.) for p in range(n)]

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

    ### Work through the parameterisations, calculating their contributions
    ### to derivatives d[pv]/dp and d[DeltaPsi]/dp

      # Calculate derivatives of pv wrt each parameter of the detector
      # parameterisations. All derivatives of DeltaPsi are zero for detector
      # parameters
      if self._detector_parameterisations:
        self._detector_derivatives(reflections, isel, dpv_dp,
                                   D, pv, det_param_id, exp.detector)

      # Calc derivatives of pv and DeltaPsi wrt each parameter of each beam
      # parameterisation that is present.
      if self._beam_parameterisations:
        self._beam_derivatives(reflections, isel, dpv_dp, dDeltaPsi_dp,
                               s0u, q0, e1, a, b, q1, r, D, beam_param_id)

      # Calc derivatives of pv and phi wrt each parameter of each crystal
      # orientation parameterisation that is present.
      if self._xl_orientation_parameterisations:
        self._xl_orientation_derivatives(reflections, isel, dpv_dp, dDeltaPsi_dp,
                                         h, B, D, xl_ori_param_id)

      # Now derivatives of pv and phi wrt each parameter of each crystal unit
      # cell parameterisation that is present.
      if self._xl_unit_cell_parameterisations:
        self._xl_unit_cell_derivatives(reflections, isel, dpv_dp, dDeltaPsi_dp,
                                       h, U, D, xl_uc_param_id)

      # calculate positional derivatives from d[pv]/dp
      dX_dp, dY_dp = self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, dpv_dp)

    return (dX_dp, dY_dp, dDeltaPsi_dp)


  def _detector_derivatives(self, reflections, isel, dpv_dp,
                            D, pv, det_param_id, detector):
    """helper function to extend the derivatives lists by derivatives of the
    detector parameterisations"""

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

  def _beam_derivatives(self, reflections, isel, dpv_dp, dDeltaPsi_dp,
                        s0u, q0, e1, a, b, q1, r, D, beam_param_id):
    """helper function to extend the derivatives lists by derivatives of the
    beam parameterisations"""

    # loop over all the beam parameterisations, even though we are only setting
    # values for one of them. We still need to move the _iparam pointer for the
    # others.
    for ibp, bp in enumerate(self._beam_parameterisations):

      # Calculate gradients only for the correct beam parameterisation
      if ibp == beam_param_id:

        # get the derivatives of the beam vector wrt the parameters
        ds0_dbeam_p = bp.get_ds_dp()

        # select indices for the experiment of interest
        sub_s0u = s0u.select(isel)
        sub_q0 = q0.select(isel)
        sub_e1 = e1.select(isel)
        sub_a = a.select(isel)
        sub_b = b.select(isel)
        sub_q1 = q1.select(isel)
        sub_r = r.select(isel)
        sub_D = D.select(isel)

        # loop through the parameters
        for der in ds0_dbeam_p:

          # calculate the derivative of DeltaPsi for this parameter
          de1 = sub_q0.cross(der)
          ds0 = flex.mat3_double(len(sub_e1), der.elems)
          dc0 = ds0.cross(sub_e1) + sub_s0u.cross(de1)
          dr = sub_b * dc0 - sub_a * ds0
          dq1 = sub_q0.cross(de1)

          y = sub_r.dot(sub_q1)
          x = sub_r.dot(sub_q0)
          dy = dr.dot(sub_q1) + sub_r.dot(dq1)
          dx = dr.dot(sub_q0)

          dDelPsi = (y * dx - x * dy) / (x**2 + y**2)

          # calculate the derivative of pv for this parameter
          dpv = sub_D * der

          # set values in the correct gradient arrays
          dDeltaPsi_dp[self._iparam].set_selected(isel, dDelPsi)
          dpv_dp[self._iparam].set_selected(isel, dpv)

          # increment the parameter index pointer
          self._iparam += 1

      # For any other beam parameterisations, leave derivatives as zero
      else:

        # just increment the pointer
        self._iparam += bp.num_free()

    return

  def _xl_orientation_derivatives(self, reflections, isel, dpv_dp, dDeltaPsi_dp,
                                  h, B, D, xl_ori_param_id):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations"""

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
        sub_B = B.select(isel)
        sub_D = D.select(isel)

        # loop through the parameters
        for der in dU_dxlo_p:

          der_mat = flex.mat3_double(len(sub_B), der.elems)
          # calculate the derivative of r for this parameter
          dq = der_mat * sub_B * sub_h

          # calculate the derivative of DeltaPsi for this parameter
          # FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
          dDelPsi = 0.0

          # calculate the derivative of pv for this parameter
          dpv = sub_D * dq

          # set values in the correct gradient arrays
          dDeltaPsi_dp[self._iparam].set_selected(isel, dDelPsi)
          dpv_dp[self._iparam].set_selected(isel, dpv)

          # increment the parameter index pointer
          self._iparam += 1

      # For any other xl orientation parameterisations, leave derivatives as zero
      else:

        # just increment the pointer
        self._iparam += xlop.num_free()

    return

  def _xl_unit_cell_derivatives(self, reflections, isel, dpv_dp, dDeltaPsi_dp,
                                h, U, D, xl_uc_param_id):
    """helper function to extend the derivatives lists by
    derivatives of the crystal unit cell parameterisations"""


    for ixlucp, xlucp in enumerate(self._xl_unit_cell_parameterisations):

      # Calculate gradients only for the correct xl unit cell parameterisation
      if ixlucp == xl_uc_param_id:

        # get derivatives of the B matrix wrt the parameters
        dB_dxluc_p = xlucp.get_ds_dp()

        # select indices for the experiment of interest
        sub_h = h.select(isel)
        sub_U = U.select(isel)
        sub_D = D.select(isel)

        # loop through the parameters
        for der in dB_dxluc_p:

          der_mat = flex.mat3_double(len(sub_U), der.elems)
          # calculate the derivative of r for this parameter
          dq = sub_U * der_mat * sub_h

          # calculate the derivative of DeltaPsi for this parameter
          # FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
          dDelPsi = 0.0

          # calculate the derivative of pv for this parameter
          dpv = sub_D * dq

          # set values in the correct gradient arrays
          dDeltaPsi_dp[self._iparam].set_selected(isel, dDelPsi)
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

    # precalculate for efficiency
    w_inv = 1/w
    u_w_inv = u * w_inv
    v_w_inv = v * w_inv

    dX_dp = []
    dY_dp = []

    for der in dpv_dp:
      du_dp, dv_dp, dw_dp = der.parts()

      dX_dp.append(w_inv * (du_dp - dw_dp * u_w_inv))
      dY_dp.append(w_inv * (dv_dp - dw_dp * v_w_inv))

    return dX_dp, dY_dp

