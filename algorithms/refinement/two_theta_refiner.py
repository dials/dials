#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Versions of refinement classes for two theta refinement of the unit cell"""

from __future__ import division
from dials.array_family import flex
from scitbx.math import angle_derivative_wrt_vectors
from math import sqrt, pi
from logging import info

from dials.algorithms.refinement.reflection_manager import ReflectionManager
from dials.algorithms.refinement.target import Target
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
      PredictionParameterisation

# constants
RAD2DEG = 180. / pi
DEG2RAD = pi / 180.

class ConstantTwoThetaWeightingStrategy(object):

  def calculate_weights(self, reflections):

    reflections['2theta.weights'] = flex.double(len(reflections), 1)
    return reflections

def calc_2theta(reflections, experiments):
  """Calculate and return 2theta angles in radians"""

  twotheta = flex.double(len(reflections), 0.)
  for iexp, exp in enumerate(experiments):
    isel = reflections['id'] == iexp
    s0 = exp.beam.get_s0()
    s1 = reflections['s1'].select(isel)
    twotheta.set_selected(isel, s1.angle(s0))
  return twotheta

class TwoThetaReflectionManager(ReflectionManager):
  _weighting_strategy = ConstantTwoThetaWeightingStrategy()

  def __init__(self, *args, **kwargs):

    # call base __init__
    super(TwoThetaReflectionManager, self).__init__(*args, **kwargs)

    # set observed 2theta angles
    self._reflections['2theta_obs.rad'] = calc_2theta(self._reflections,
      self._experiments)

    return

  def print_stats_on_matches(self):

    l = self.get_matches()
    nref = len(l)

    from libtbx.table_utils import simple_table
    from scitbx.math import five_number_summary
    twotheta_resid = l['2theta_resid']
    w_2theta = l['2theta.weights']

    msg = "\nSummary statistics for {0} observations".format(nref) +\
          " matched to predictions:"
    header = ["", "Min", "Q1", "Med", "Q3", "Max"]
    rows = []
    try:
      row_data = five_number_summary(twotheta_resid)
      rows.append(["2theta_c - 2theta_o (deg)"] + ["%.4g" % (e * RAD2DEG) for e in row_data])
      row_data = five_number_summary(w_2theta)
      rows.append(["2theta weights"] + ["%.4g" % (e * DEG2RAD**2) for e in row_data])
      st = simple_table(rows, header)
    except IndexError:
      # zero length reflection list
      warning("Unable to calculate summary statistics for zero observations")
      return
    info(msg)
    info(st.format())
    info("")

class TwoThetaTarget(Target):
  _grad_names = ['d2theta_dp']
  rmsd_names = ["RMSD_2theta"]
  rmsd_units = ["rad"]

  def __init__(self, experiments, reflection_predictor, ref_man,
               prediction_parameterisation):
    Target.__init__(self, experiments, reflection_predictor, ref_man,
                    prediction_parameterisation)

    # set the single cutoff for 2theta residual to essentially zero
    self._binsize_cutoffs = [1.e-6]

    # predict reflections and finalise reflection manager
    self.predict()
    self._reflection_manager.finalise()

    return

  def predict(self):
    """perform reflection prediction for the working reflections and update the
    reflection manager"""

    # get the matches
    reflections = self._reflection_manager.get_obs()

    # reset the 'use' flag for all observations
    self._reflection_manager.reset_accepted_reflections()

    # standard centroid prediction
    reflections = self._predict_core(reflections)

    # calculate two theta values and residuals
    twotheta = calc_2theta(reflections, self._experiments)
    reflections['2theta_resid'] = twotheta - reflections['2theta_obs.rad']
    reflections['2theta_resid2'] = reflections['2theta_resid']**2

    # set used_in_refinement flag to all those that had predictions
    mask = reflections.get_flags(reflections.flags.predicted)
    reflections.set_flags(mask, reflections.flags.used_in_refinement)

    # collect the matches
    self.update_matches(force=True)

    return

  # XXXX FIXME In the base class this method is hardcoded to expect three
  # types of residual: in X, Y, and Z. Better to make the base class behave
  # with any number of types and remove this method
  @staticmethod
  def _build_jacobian(d2theta_dp, nelem=None, nparam=None):
    """construct Jacobian from lists of gradient vectors. This method may be
    overridden for the case where these vectors use sparse storage"""

    jacobian = flex.double(flex.grid(nelem, nparam))
    # loop over parameters
    for i in range(nparam):
      col = d2theta_dp[i]
      jacobian.matrix_paste_column_in_place(col, i)

    return jacobian

  @staticmethod
  def _extract_residuals_and_weights(matches):

    # return residuals and weights as 1d flex.double vectors
    residuals = matches['2theta_resid']

    weights = matches['2theta.weights']

    return residuals, weights

  @staticmethod
  def _extract_squared_residuals(matches):

    residuals2 = matches['2theta_resid2']

    return residuals2

  def _rmsds_core(self, reflections):
    """calculate unweighted RMSDs for the specified reflections"""

    resid_2theta = flex.sum(reflections['2theta_resid2'])
    n = len(reflections)

    rmsds = (sqrt(resid_2theta / n), )
    return rmsds

  def achieved(self):
    """RMSD criterion for target achieved """
    r = self._rmsds if self._rmsds else self.rmsds()

    # reset cached rmsds to avoid getting out of step
    self._rmsds = None

    if (r[0] < self._binsize_cutoffs[0]):
      return True
    return False

class TwoThetaPredictionParameterisation(PredictionParameterisation):
  _grad_names = ("d2theta_dp",)

  def __init__(self, *args, **kwargs):
    super(TwoThetaPredictionParameterisation, self).__init__(*args, **kwargs)
    # check that only the unit cell is parameterised
    assert not self._detector_parameterisations
    assert not self._beam_parameterisations
    assert not self._xl_orientation_parameterisations
    return

  def _local_setup(self, reflections):
    # r is the reciprocal lattice vector, in the lab frame
    self._phi_calc = reflections['xyzcal.mm'].parts()[2]
    q = self._fixed_rotation * (self._UB * self._h)
    self._r = q.rotate_around_origin(self._axis, self._phi_calc)

    # All of the derivatives of phi have a common denominator, given by
    # (e X r).s0, where e is the rotation axis. Calculate this once, here.
    self._e_X_r = self._axis.cross(self._r)
    self._e_r_s0 = (self._e_X_r).dot(self._s0)

    # catch small values of e_r_s0
    e_r_s0_mag = flex.abs(self._e_r_s0)
    try:
      assert flex.min(e_r_s0_mag) > 1.e-6
    except AssertionError as e:
      imin = flex.min_index(e_r_s0_mag)
      print "(e X r).s0 too small:"
      print "for", (e_r_s0_mag <= 1.e-6).count(True), "reflections"
      print "out of", len(e_r_s0_mag), "total"
      print "such as", reflections['miller_index'][imin]
      print "with scattering vector", reflections['s1'][imin]
      print "where r =", self._r[imin]
      print "e =", self._axis[imin]
      print "s0 =", self._s0[imin]
      print ("this reflection forms angle with the equatorial plane "
             "normal:")
      vecn = matrix.col(self._s0[imin]).cross(matrix.col(self._axis[imin])).normalize()
      print matrix.col(reflections['s1'][imin]).accute_angle(vecn)
      raise e

    # we want the wavelength
    self._wavelength = 1. / self._s0.norms()

    return

  def _xl_unit_cell_derivatives(self, isel, parameterisation=None,
    reflections=None):

    # Get required data
    h = self._h.select(isel)
    B = self._B.select(isel)
    wl = self._wavelength.select(isel)

    # get derivatives of the B matrix wrt the parameters
    dB_dxluc_p = [None if der is None else flex.mat3_double(len(isel), der.elems) \
                  for der in parameterisation.get_ds_dp(use_none_as_null=True)]

    d2theta_dp = []

    # loop through the parameters
    for der in dB_dxluc_p:

      if der is None:
        d2theta_dp.append(None)
        continue

      r0 = B * h
      dr0 = der * h
      r0len = r0.norms()
      dr0len = dr0.dot(r0) / r0len

      # 2theta = 2 * arcsin( |r0| / (2 * |s0| ) )
      sintheta = 0.5 * r0len * wl
      fac = 1.0 / flex.sqrt(flex.double(len(wl), 1.0) - sintheta**2)
      val = fac * wl * dr0len

      d2theta_dp.append(val)

    return d2theta_dp

  def _grads_xl_unit_cell_loop(self, reflections, results, callback=None):
    """Loop over all crystal unit cell parameterisations, calculate gradients
    and extend the results"""

    # loop over the crystal unit cell parameterisations
    for xlucp in self._xl_unit_cell_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in xlucp.get_experiment_ids():
        isel.extend(self._experiment_to_idx[exp_id])

      # Extend derivative vectors for this crystal unit cell parameterisation
      results = self._extend_gradient_vectors(results, self._nref, xlucp.num_free(),
        keys=self._grad_names)

      if len(isel) == 0:
        # if no reflections are in this experiment, skip calculation of
        # gradients, but must still process null gradients by a callback
        if callback is not None:
          for iparam in xrange(xlucp.num_free()):
            results[self._iparam] = callback(results[self._iparam])
            self._iparam += 1
        else:
          self._iparam += xlucp.num_free()
        continue

      d2theta_dp =  self._xl_unit_cell_derivatives(isel,
        parameterisation=xlucp, reflections=reflections)

      for d2theta in d2theta_dp:
        if d2theta is not None:
          results[self._iparam][self._grad_names[0]].set_selected(isel, d2theta)

        # increment the parameter index pointer
        self._iparam += 1

    return results
