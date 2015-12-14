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

# dials imports
from dials.algorithms.refinement.target import Target, SparseGradientsMixin

# constants
TWO_PI = 2.0 * pi

class LeastSquaresStillsResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y and minimum rotation to the Ewald
  sphere, DeltaPsi. Terminates refinement on achieved rmsd (or on intrisic
  convergence of the chosen minimiser)"""

  _grad_names = ['dX_dp', 'dY_dp', 'dDeltaPsi_dp']
  rmsd_names = ["RMSD_X", "RMSD_Y", "RMSD_DeltaPsi"]
  rmsd_units = ["mm", "mm", "rad"]

  def __init__(self, experiments, reflection_predictor, ref_man,
               prediction_parameterisation, restraints_parameterisation,
               frac_binsize_cutoff=0.33333,
               absolute_cutoffs=None,
               gradient_calculation_blocksize=None):

    Target.__init__(self, experiments, reflection_predictor, ref_man,
                    prediction_parameterisation, gradient_calculation_blocksize)

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

  def  _predict_core(self, reflections):

    """perform prediction for the specified reflections"""
    # update the reflection_predictor with the scan-independent part of the
    # current geometry
    self._reflection_predictor.update()

    # do prediction (updates reflection table in situ).
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

    return reflections

  def predict_for_reflection_table(self, reflections):
    """perform prediction for all reflections in the supplied table"""

    # predict
    return self._predict_core(reflections)

  @staticmethod
  def _extract_residuals_and_weights(matches):

    # return residuals and weights as 1d flex.double vectors
    residuals = flex.double.concatenate(matches['x_resid'],
                                        matches['y_resid'])
    residuals.extend(matches['delpsical.rad'])

    weights, w_y, _ = matches['xyzobs.mm.weights'].parts()
    w_delpsi = matches['delpsical.weights']
    weights.extend(w_y)
    weights.extend(w_delpsi)

    return residuals, weights

  @staticmethod
  def _extract_squared_residuals(matches):

    residuals2 = flex.double.concatenate(matches['x_resid2'],
                                         matches['y_resid2'])
    residuals2.extend(matches['delpsical2'])

    return residuals2

  def _rmsds_core(self, reflections):
    """calculate unweighted RMSDs"""

    resid_x = flex.sum(reflections['x_resid2'])
    resid_y = flex.sum(reflections['y_resid2'])
    resid_z = flex.sum(reflections['delpsical2'])
    n = len(reflections)
    rmsds = (sqrt(resid_x / n),
             sqrt(resid_y / n),
             sqrt(resid_z / n))
    return rmsds

  def achieved(self):
    """RMSD criterion for target achieved """
    r = self._rmsds if self._rmsds else self.rmsds()

    # reset cached rmsds to avoid getting out of step
    self._rmsds = None

    # only use RMSD_X and RMSD_Y
    if (r[0] < self._binsize_cutoffs[0] and
        r[1] < self._binsize_cutoffs[1]):
      return True
    return False

class LeastSquaresStillsResidualWithRmsdCutoffSparse(
    SparseGradientsMixin, LeastSquaresStillsResidualWithRmsdCutoff):
  """A version of the LeastSquaresStillsResidualWithRmsdCutoff Target that
  uses a sparse matrix data structure for memory efficiency when there are a
  large number of Experiments"""

  pass
