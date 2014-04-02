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
#from dials_refinement_helpers_ext import *
from dials.algorithms.refinement.parameterisation.prediction_parameters_new import \
    XYPhiPredictionParameterisation

from dials.array_family import flex

class VaryingCrystalPredictionParameterisation(XYPhiPredictionParameterisation):
  """Support crystal parameterisations that vary with time (via the proxy of
  "observed image number")"""

  def compose(self, reflections):
    """Compose scan-varying crystal parameterisations at the specified image
    number, for the specified experiment, for all reflections. Put the U, B and
    UB matrices in the reflection table, and cache the derivatives."""

    nref = len(reflections)
    # set columns if needed
    if not reflections.has_key('u_matrix'):
      reflections['u_matrix'] = flex.mat3_double(nref)
    if not reflections.has_key('b_matrix'):
      reflections['b_matrix'] = flex.mat3_double(nref)

    # set up arrays to store derivatives
    num_free_U_params = sum([e.num_free() for e in self._xl_orientation_parameterisations])
    num_free_B_params = sum([e.num_free() for e in self._xl_unit_cell_parameterisations])
    self._dU_dp = [flex.mat3_double(nref) for i in range(num_free_U_params)]
    self._dB_dp = [flex.mat3_double(nref) for i in range(num_free_B_params)]

    for iexp, exp in enumerate(self._experiments):

      # select the reflections of interest
      sel = reflections['id'] == iexp
      isel = sel.iselection()

      # get their frame numbers
      obs_image_numbers = (reflections['xyzobs.px.value'].parts()[2]).select(isel)

      # identify which crystal parameterisations to use for this experiment
      param_set = self._exp_to_param[iexp]
      xl_ori_param_id = param_set.xl_ori_param
      xl_uc_param_id = param_set.xl_uc_param
      xl_op = self._xl_orientation_parameterisations[param_set.xl_ori_param]
      xl_ucp = self._xl_unit_cell_parameterisations[param_set.xl_uc_param]

      # get state and derivatives for each reflection
      for i, frame in zip(isel, obs_image_numbers):

        # compose the models
        xl_op.compose(frame)
        xl_ucp.compose(frame)

        # set states
        row = {'u_matrix':xl_op.get_state().elems,
               'b_matrix':xl_ucp.get_state().elems}
        reflections[i] = row

        # set derivatives of the states
        for j, dU in enumerate(xl_op.get_ds_dp()):
          self._dU_dp[j][i] = dU

        for j, dB in enumerate(xl_ucp.get_ds_dp()):
          self._dB_dp[j][i] = dB

    # set the UB matrices for prediction
    reflections['ub_matrix'] = reflections['u_matrix'] * reflections['b_matrix']

    return

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
    n = len(reflections)
    D = flex.mat3_double(n)
    s0 = flex.vec3_double(n)
    #U = flex.mat3_double(n)
    #B = flex.mat3_double(n)
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

      # axis array
      if exp.goniometer:
        axis.set_selected(isel, exp.goniometer.get_rotation_axis())

    # for the scan varying parameterisation, U and B are set in the reflections
    # table by the compose method
    U = reflections['u_matrix']
    B = reflections['b_matrix']

    return self._get_gradients_core(reflections, D, s0, U, B, axis)

  # overloaded for the scan-varying case
  def _xl_orientation_derivatives(self, reflections, isel, dpv_dp, dphi_dp, axis, phi_calc, h, s1, \
                                         e_X_r, e_r_s0, B, D, xl_ori_param_id):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations."""

    # loop over all the crystal orientation parameterisations, even though we
    # are only setting values for one of them. We still need to move the _iparam
    # pointer for the others.
    for ixlop, xlop in enumerate(self._xl_orientation_parameterisations):

      # Calculate gradients only for the correct xl orientation parameterisation
      if ixlop == xl_ori_param_id:

        # get derivatives of the U matrix wrt the parameters
        #dU_dxlo_p = self._exp_to_xl_derivatives[iexp].dU_dp
        dU_dxlo_p = [e.select(isel) for e in self._dU_dp]

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

          assert len(der) == len(isel)
          # each der is a flex.mat3_double array
          tmp = der * sub_B * sub_h
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

  # overloaded for the scan-varying case
  def _xl_unit_cell_derivatives(self, reflections, isel, dpv_dp, dphi_dp, axis, phi_calc, h, s1, \
                                         e_X_r, e_r_s0, U, D, xl_uc_param_id):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations."""

    for ixlucp, xlucp in enumerate(self._xl_unit_cell_parameterisations):

      # Calculate gradients only for the correct xl unit cell parameterisation
      if ixlucp == xl_uc_param_id:

        # get derivatives of the B matrix wrt the parameters
        #dB_dxluc_p = self._exp_to_xl_derivatives[iexp].dB_dp
        dB_dxluc_p = [e.select(isel) for e in self._dB_dp]

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

          assert len(der) == len(isel)
          # each der is a flex.mat3_double array
          # calculate the derivative of r for this parameter
          tmp = sub_U * der * sub_h
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
