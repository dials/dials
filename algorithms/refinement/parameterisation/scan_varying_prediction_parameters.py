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
from dials_refinement_helpers_ext import *
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation

# A helper bucket to store cached values
from collections import namedtuple
ModelCache = namedtuple('ModelCache', ['D_mats', 's0', 'axis'])

class VaryingCrystalPredictionParameterisation(XYPhiPredictionParameterisation):
  """Support crystal parameterisations that vary with time (via the proxy of
  "observed image number")"""

  _obs_image_number = None

  def prepare(self, panel_id=0):
    """Cache required quantities that are not dependent on hkl"""

    # Same as prepare for the parent class except we don't get
    # U and B from the model
    self._cache = []
    for e in self._experiments:

      D_mats=[matrix.sqr(p.get_D_matrix()) for p in e.detector]
      s0 = matrix.col(e.beam.get_s0())
      if e.goniometer:
        axis = matrix.col(e.goniometer.get_rotation_axis())
      else:
        axis = None

      self._cache.append(ModelCache(D_mats, s0, axis=axis))

    return

  def compose(self, obs_image_number):
    """Compose scan-varying crystal parameterisations at the specified
    image number"""

    self._obs_image_number = obs_image_number
    xl_op = self._xl_orientation_parameterisations[0]
    xl_ucp = self._xl_unit_cell_parameterisations[0]
    xl_op.compose(obs_image_number)
    xl_ucp.compose(obs_image_number)

  def get_UB(self, obs_image_number):
    """Extract the setting matrix from the contained scan
    dependent crystal parameterisations at specified image number"""

    if obs_image_number != self._obs_image_number:
      self.compose(obs_image_number)

    UB = self._xl_orientation_parameterisations[0].get_state() * \
         self._xl_unit_cell_parameterisations[0].get_state()
    return UB

  def get_gradients(self, h, s, phi, panel_id, obs_image_number,
                    experiment_id=0):

    #self.prepare()
    if obs_image_number != self._obs_image_number:
      self.compose(obs_image_number)

    # extract the right models for the requested experiment
    self._D = self._cache[experiment_id].D_mats[panel_id]
    self._s0 = self._cache[experiment_id].s0
    self._axis = self._cache[experiment_id].axis

    return self._get_gradients_core(h, s, phi, panel_id, experiment_id)

  def _get_gradients_core(self, h, s, phi, panel_id, experiment_id):

    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    that reflects at rotation angle phi with scattering vector s
    that intersects panel panel_id. That is, calculate dX/dp, dY/dp
    and dphi/dp. Scan-varying parameters (for the crystal) are
    evaluated at obs_image_number"""

    # NB prepare and compose must be called first

    ### Calculate various quantities of interest for this reflection

    R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)

    # Get U and B at the composed image number obs_image_number. Assume
    # there is only one parameterisation of each type
    xl_op = self._xl_orientation_parameterisations[0]
    xl_ucp = self._xl_unit_cell_parameterisations[0]
    U = xl_op.get_state()
    B = xl_ucp.get_state()

    # pv is the 'projection vector' for the reflection s.
    s = matrix.col(s)
    pv = self._D * s
    # r is the reciprocal lattice vector, in the lab frame
    r = R * U * B * h

    # All of the derivatives of phi have a common denominator, given by
    # (e X r).s0, where e is the rotation axis. Calculate this once, here.
    e_X_r = self._axis.cross(r)
    e_r_s0 = (e_X_r).dot(self._s0)

    # Note that e_r_s0 -> 0 when the rotation axis, beam vector and
    # relp are coplanar. This occurs when a reflection just touches
    # the Ewald sphere.
    #
    # There is a relationship between e_r_s0 and zeta_factor.
    # Uncommenting the code below shows that
    # s0.(e X r) = zeta * |s X s0|

    #from dials.algorithms.reflection_basis import zeta_factor
    #from libtbx.test_utils import approx_equal
    #z = zeta_factor(self._axis, self._s0, s)
    #ss0 = (s.cross(self._s0)).length()
    #assert approx_equal(e_r_s0, z * ss0)

    # catch small values of e_r_s0
    try:
      assert abs(e_r_s0) > 1.e-6
    except AssertionError as e:
      print "(e X r).s0 too small:", e_r_s0
      print "for reflection", h
      print "with scattering vector", s
      print "where r =", r
      print "e =",matrix.col(self._axis)
      print "s0 =",self._s0
      print "U(t) =",U
      print ("this reflection forms angle with the equatorial plane "
             "normal:")
      vecn = self._s0.cross(self._axis).normalize()
      print s.accute_angle(vecn)
      raise e

    # identify which parameterisations to use
    param_set = self._exp_to_param[experiment_id]
    beam_param_id = param_set.beam_param
    xl_ori_param_id = param_set.xl_ori_param
    xl_uc_param_id = param_set.xl_uc_param
    det_param_id = param_set.det_param

    ### Work through the parameterisations, calculating their contributions
    ### to derivatives d[pv]/dp and d[phi]/dp

    # Set up the lists of derivatives
    dpv_dp = []
    dphi_dp = []

    # Calculate derivatives of pv wrt each parameter of the detector
    # parameterisations. All derivatives of phi are zero for detector
    # parameters
    if self._detector_parameterisations:
      self._detector_derivatives(dpv_dp, dphi_dp, pv, panel_id, det_param_id)

    # Calc derivatives of pv and phi wrt each parameter of each beam
    # parameterisation that is present.
    if self._beam_parameterisations:
      self._beam_derivatives(dpv_dp, dphi_dp, r, e_X_r, e_r_s0, beam_param_id)

    # Calc derivatives of pv and phi wrt each parameter of each
    # scan-varying crystal orientation parameterisation
    if self._xl_orientation_parameterisations:
      self._xl_orientation_derivatives(dpv_dp, dphi_dp, \
              self._obs_image_number, B, R, h, s, e_X_r, e_r_s0, xl_ori_param_id)

    # Now derivatives of pv and phi wrt each parameter of each
    # scan-varying crystal unit cell parameterisation
    if self._xl_unit_cell_parameterisations:
      self._xl_unit_cell_derivatives(dpv_dp, dphi_dp, \
              self._obs_image_number, U, R, h, s, e_X_r, e_r_s0, xl_uc_param_id)

    # calculate positional derivatives from d[pv]/dp
    pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e)
                for e in dpv_dp]
    dX_dp, dY_dp = zip(*pos_grad)

    return zip(dX_dp, dY_dp, dphi_dp)

  def _xl_orientation_derivatives(self, dpv_dp, dphi_dp, \
          obs_image_number, B, R, h, s, e_X_r, e_r_s0, xl_ori_param_id):
    """Adds calculation at obs_image_number for scan-varying
    parameters"""

    for ixlo, xlo in enumerate(self._xl_orientation_parameterisations):

      # Calculate gradients only for the correct xl orientation parameterisation
      if ixlo == xl_ori_param_id:
        dU_dxlo_p = xlo.get_ds_dp(obs_image_number)

        dr_dxlo_p = [R * dU_dxlo_p[i] * B * h
                     for i in range(len(dU_dxlo_p))]

        dphi_dxlo_p = [- der.dot(s) / e_r_s0 for der in dr_dxlo_p]

        dpv_dxlo_p = [self._D * (dr_dxlo_p[i] + e_X_r * dphi_dxlo_p[i])
                      for i in range(len(dphi_dxlo_p))]

      # For any other xl orientation parameterisations, set derivatives to zero
      else:
        dphi_dxlo_p = [0.] * len(xlo.num_free())
        dpv_dxlo_p = [matrix.col((0., 0., 0.))] * len(xlo.num_free())

      dpv_dp.extend(dpv_dxlo_p)
      dphi_dp.extend(dphi_dxlo_p)

    return

  def _xl_unit_cell_derivatives(self, dpv_dp, dphi_dp, \
          obs_image_number, U, R, h, s, e_X_r, e_r_s0, xl_uc_param_id):
    """Adds calculation at obs_image_number for scan-varying
    parameters"""

    for ixluc, xluc in enumerate(self._xl_unit_cell_parameterisations):

      # Calculate gradients only for the correct xl unit cell parameterisation
      if ixluc == xl_uc_param_id:
        dB_dxluc_p = xluc.get_ds_dp(obs_image_number)

        dr_dxluc_p = [R * U * dB_dxluc_p[i] * h for i
                          in range(len(dB_dxluc_p))]

        dphi_dxluc_p = [- der.dot(s) / e_r_s0 for der in dr_dxluc_p]

        dpv_dxluc_p = [self._D * (dr_dxluc_p[i] + e_X_r * dphi_dxluc_p[i])
                       for i in range(len(dr_dxluc_p))]

      # For any other xl unit cell parameterisations, set derivatives to zero
      else:
        dphi_dxluc_p = [0.] * len(xlo.num_free())
        dpv_dxluc_p = [matrix.col((0., 0., 0.))] * len(xlo.num_free())

      dpv_dp.extend(dpv_dxluc_p)
      dphi_dp.extend(dphi_dxluc_p)

    return

  def _calc_dX_dp_and_dY_dp_from_dpv_dp(self, pv, der):
    """helper function to calculate positional derivatives from dpv_dp using
    the quotient rule"""
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
