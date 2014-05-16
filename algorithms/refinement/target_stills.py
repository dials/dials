#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# python and cctbx imports
from __future__ import division
from math import pi, sqrt
from cctbx.array_family import flex
import random

# dials imports
from dials.algorithms.refinement.target import Target

# constants
TWO_PI = 2.0 * pi

class LeastSquaresStillsResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y and minimum rotation to the Ewald
  sphere, DeltaPsi. Terminates refinement on achieved rmsd (or on intrisic
  convergence of the chosen minimiser)"""

  rmsd_names = ["RMSD_X", "RMSD_Y"]

  def __init__(self, experiments, reflection_predictor, ref_man,
               prediction_parameterisation,
               frac_binsize_cutoff=0.33333,
               absolute_cutoffs=None,
               jacobian_max_nref=None):

    Target.__init__(self, experiments, reflection_predictor, ref_man,
                    prediction_parameterisation, jacobian_max_nref)

    # Set up the RMSD achieved criterion. For simplicity, we take models from
    # the first Experiment only. If this is not appropriate for refinement over
    # all experiments then absolute cutoffs should be used instead.
    detector = experiments[0].detector
    if not absolute_cutoffs:
      pixel_sizes = [p.get_pixel_size() for p in detector]
      min_px_size_x = min(e[0] for e in pixel_sizes)
      min_px_size_y = min(e[1] for e in pixel_sizes)
      self._binsize_cutoffs = [min_px_size_x * frac_binsize_cutoff,
                               min_px_size_y * frac_binsize_cutoff]
    else:
      self._binsize_cutoffs = absolute_cutoffs[:2]

    # predict reflections and finalise reflection manager
    self.predict()
    self._reflection_manager.finalise()

    return


  def predict(self):
    """perform reflection prediction and update the reflection manager"""

    # update the reflection_predictor with the scan-independent part of the
    # current geometry
    self._reflection_predictor.update()

    # reset the 'use' flag for all observations
    self._reflection_manager.reset_accepted_reflections()

    # do prediction (updates reflection table in situ).
    reflections = self._reflection_manager.get_obs()
    self._reflection_predictor.predict(reflections)

    x_obs, y_obs, _ = reflections['xyzobs.mm.value'].parts()
    delpsi = reflections['delpsical.rad']
    x_calc, y_calc, _ = reflections['xyzcal.mm'].parts()

    # calculate residuals and assign columns
    reflections['x_resid'] = x_calc - x_obs
    reflections['x_resid2'] = reflections['x_resid']**2
    reflections['y_resid'] = y_calc - y_obs
    reflections['y_resid2'] = reflections['y_resid']**2
    reflections['delpsical2'] = reflections['delpsical.rad']**2

    # set used_in_refinement flag to all those that had predictions
    mask = reflections.get_flags(reflections.flags.predicted)
    reflections.set_flags(mask, reflections.flags.used_in_refinement)

    # collect the matches
    self._matches = self._reflection_manager.get_matches()

    return

  def predict_for_reflection_table(self, reflections):
    """perform prediction for all reflections in the supplied table"""

    self._reflection_predictor.update()
    self._reflection_predictor.predict(reflections)
    return reflections

  def compute_residuals_and_gradients(self):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    dX_dp, dY_dp, dDPsi_dp = self.calculate_gradients()

    # return residuals and weights as 1d flex.double vectors
    nelem = len(self._matches) * 3
    nparam = len(self._prediction_parameterisation)
    residuals = flex.double.concatenate(self._matches['x_resid'],
                                        self._matches['y_resid'])
    residuals.extend(self._matches['delpsical.rad'])
    #jacobian_t = flex.double(flex.grid(
    #    len(self._prediction_parameterisation), nelem))
    weights, w_y, _ = self._matches['xyzobs.mm.weights'].parts()
    w_delpsi = self._matches['delpsical.weights']
    weights.extend(w_y)
    weights.extend(w_delpsi)

    jacobian = flex.double(flex.grid(nelem, nparam))
    # loop over parameters
    for i in range(nparam):
      dX, dY, dDPsi = dX_dp[i], dY_dp[i], dDPsi_dp[i]
      col = flex.double.concatenate(dX, dY)
      col.extend(dDPsi)
      jacobian.matrix_paste_column_in_place(col, i)

    return(residuals, jacobian, weights)

  def compute_functional_and_gradients(self):
    """calculate the value of the target function and its gradients"""

    self._nref = len(self._matches)

    # This is a hack for the case where nref=0. This should not be necessary
    # if bounds are provided for parameters to stop the algorithm exploring
    # unreasonable regions of parameter space where no predictions exist.
    # Unfortunately the L-BFGS line search does make such extreme trials.
    if self._nref == 0:
      return 1.e12, [1.] * len(self._prediction_parameterisation)

    # extract columns from the table
    x_resid = self._matches['x_resid']
    y_resid = self._matches['y_resid']
    delpsi = self._matches['delpsical.rad']
    x_resid2 = self._matches['x_resid2']
    y_resid2 = self._matches['y_resid2']
    delpsical2 = self._matches['delpsical2']
    w_x, w_y, _ = self._matches['xyzobs.mm.weights'].parts()
    w_delpsi = self._matches['delpsical.weights']

    # calculate target function
    temp = w_x * x_resid2 + w_y * y_resid2 + w_delpsi * delpsical2
    L = 0.5 * flex.sum(temp)

    # prepare list of gradients
    dL_dp = [0.] * len(self._prediction_parameterisation)

    dX_dp, dY_dp, dDPsi_dp = self.calculate_gradients()

    w_x_x_resid = w_x * x_resid
    w_y_y_resid = w_y * y_resid
    w_delpsi_delpsi = w_delpsi * delpsi

    for i in range(len(self._prediction_parameterisation)):
      dX, dY, dDPsi = dX_dp[i], dY_dp[i], dDPsi_dp[i]
      temp = w_x_x_resid * dX + w_y_y_resid * dY + w_delpsi_delpsi * dDPsi
      dL_dp[i] = flex.sum(temp)

    return (L, dL_dp)

  def curvatures(self):
    """First order approximation to the diagonal of the Hessian based on the
    least squares form of the target"""

    # relies on compute_functional_and_gradients being called first
    dX_dp, dY_dp, dDPsi_dp = self._gradients
    w_x, w_y, _ = self._matches['xyzobs.mm.weights'].parts()
    w_delpsi = self._matches['delpsical.weights']

    # This is a hack for the case where nref=0. This should not be necessary
    # if bounds are provided for parameters to stop the algorithm exploring
    # unreasonable regions of parameter space where there are no predictions
    if self._nref == 0:
      return [1.] * len(self._prediction_parameterisation)

    # prepare lists of gradients and curvatures
    curv = [0.] * len(self._prediction_parameterisation)

    # for each reflection, get the approximate curvatures wrt each parameter
    for i in range(len(self._prediction_parameterisation)):
      dX, dY, dDPsi = dX_dp[i], dY_dp[i], dDPsi_dp[i]
      temp = w_x * dX**2 + w_y * dY**2 + w_delpsi * dDPsi**2
      curv[i] = flex.sum(temp)

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert all([c > 0.0 for c in curv])

    return curv

  def rmsds(self):
    """calculate unweighted RMSDs"""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()

    resid_x = flex.sum(self._matches['x_resid2'])
    resid_y = flex.sum(self._matches['y_resid2'])

    # cache rmsd calculation for achieved test
    n = len(self._matches)
    self._rmsds = (sqrt(resid_x / n),
                   sqrt(resid_y / n))

    return self._rmsds

  def achieved(self):
    """RMSD criterion for target achieved """
    r = self._rmsds if self._rmsds else self.rmsds()

    # reset cached rmsds to avoid getting out of step
    self._rmsds = None

    if (r[0] < self._binsize_cutoffs[0] and
        r[1] < self._binsize_cutoffs[1]):
      return True
    return False
