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

#### DIALS imports

from dials_refinement_helpers_ext import *
from dials.algorithms.refinement.parameterisation import \
    PredictionParameterisation

class XYPredictionParameterisation(PredictionParameterisation):
  """
  Concrete class that inherits functionality of the
  PredictionParameterisation parent class and provides a detector X, Y space
  implementation of the get_gradients function.

  Untested for multiple sensor detectors.
  """

  def _get_gradients_core(self, h, s, phi, panel_id):

    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    that reflects at rotation angle phi with scattering vector s that
    intersects panel panel_id. That is, calculate dX/dp, dY/dp and
    dphi/dp"""

    ### Calculate various quantities of interest for this reflection

    # Are we dealing with rotation data or a still?
    if self._axis:
      R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)
    else:
      R = matrix.identity(3)

    # pv is the 'projection vector' for the reflection s.
    s = matrix.col(s)
    pv = self._D * s
    # r is the reciprocal lattice vector, in the lab frame
    r = R * self._UB * h

    ### Work through the parameterisations, calculating their contributions
    ### to derivatives d[pv]/dp

    # Set up the lists of derivatives
    dpv_dp = []

    # Calculate derivatives of pv wrt each parameter of the FIRST detector
    # parameterisation only.
    if self._detector_parameterisations:
      self._detector_derivatives(dpv_dp, pv, panel_id)

    # Calc derivatives of pv wrt each parameter of each beam
    # parameterisation that is present.
    if self._beam_parameterisations:
      self._beam_derivatives(dpv_dp, r)

    # Calc derivatives of pv wrt each parameter of each crystal
    # orientation parameterisation that is present.
    if self._xl_orientation_parameterisations:
      self._xl_orientation_derivatives(dpv_dp, R, h)

    # Now derivatives of pv wrt each parameter of each crystal unit
    # cell parameterisation that is present.
    if self._xl_unit_cell_parameterisations:
      self._xl_unit_cell_derivatives(dpv_dp, R, h)

    # calculate positional derivatives from d[pv]/dp
    pos_grad = [self._calc_dX_dp_and_dY_dp_from_dpv_dp(pv, e) for e in dpv_dp]
    dX_dp, dY_dp = zip(*pos_grad)

    return zip(dX_dp, dY_dp)

  def _detector_derivatives(self, dpv_dp, pv, panel_id):
    """helper function to extend the derivatives lists by
    derivatives of the detector parameterisations"""

    for idet, det in enumerate(self._detector_parameterisations):
      if idet == 0:
        dd_ddet_p = det.get_ds_dp(multi_state_elt=panel_id)
        dpv_ddet_p = [- self._D * e * pv for e in dd_ddet_p]
      else:
        dpv_ddet_p = [matrix.col((0., 0., 0.))] * len(dd_ddet_p)

      dpv_dp.extend(dpv_ddet_p)

    return

  def _beam_derivatives(self, dpv_dp, r):
    """helper function to extend the derivatives lists by
    derivatives of the beam parameterisations"""

    for src in self._beam_parameterisations:
      ds0_dsrc_p = src.get_ds_dp()
      dpv_dsrc_p = [self._D * e for e in ds0_dsrc_p]

      dpv_dp.extend(dpv_dsrc_p)

    return

  def _xl_orientation_derivatives(self, dpv_dp, R, h):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations"""

    for xlo in self._xl_orientation_parameterisations:
      dU_dxlo_p = xlo.get_ds_dp()

      dr_dxlo_p = [R * e * self._B * h for e in dU_dxlo_p]

      dpv_dxlo_p = [self._D * e for e in dr_dxlo_p]

      dpv_dp.extend(dpv_dxlo_p)

    return

  def _xl_orientation_derivatives(self, dpv_dp, R, h):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations"""

    for xlo in self._xl_orientation_parameterisations:
      dU_dxlo_p = xlo.get_ds_dp()

      dr_dxlo_p = [R * e * self._B * h for e in dU_dxlo_p]

      dpv_dxlo_p = [self._D * e for e in dr_dxlo_p]

      dpv_dp.extend(dpv_dxlo_p)

    return

  def _xl_unit_cell_derivatives(self, dpv_dp, R, h):
    """helper function to extend the derivatives lists by
    derivatives of the crystal unit cell parameterisations"""

    for xluc in self._xl_unit_cell_parameterisations:
      dB_dxluc_p = xluc.get_ds_dp()

      dr_dxluc_p = [R * self._U * e * h for e in dB_dxluc_p]

      dpv_dxluc_p = [self._D * e for e in dr_dxluc_p]

      dpv_dp.extend(dpv_dxluc_p)

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
