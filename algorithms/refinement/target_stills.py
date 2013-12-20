#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
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
from dials.algorithms.refinement.target import ReflectionManager
from dials.algorithms.refinement.target import ObsPredMatch # implicit import

# constants
TWO_PI = 2.0 * pi

class LeastSquaresXYResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y only, terminating on achieved
  rmsd (or on intrisic convergence of the chosen minimiser)"""

  rmsd_names = ["RMSD_X", "RMSD_Y"]

  def __init__(self, reflection_predictor, detector, ref_man,
               prediction_parameterisation,
               frac_binsize_cutoff=0.33333,
               absolute_cutoffs=None):

    Target.__init__(self, reflection_predictor, detector, ref_man,
                    prediction_parameterisation)

    # Set up the RMSD achieved criterion
    if not absolute_cutoffs:
      pixel_sizes = [p.get_pixel_size() for p in detector]
      min_px_size_x = min(e[0] for e in pixel_sizes)
      min_px_size_y = min(e[1] for e in pixel_sizes)
      self._binsize_cutoffs = [min_px_size_x * frac_binsize_cutoff,
                               min_px_size_y * frac_binsize_cutoff]
    else:
      self._binsize_cutoffs = absolute_cutoffs[:2]

    return

  def compute_residuals_and_gradients(self):
    """return the vector of residuals plus their gradients
    and weights for non-linear least squares methods"""

    self._matches = self._reflection_manager.get_matches()

    # return residuals and weights as 1d flex.double vectors
    nelem = len(self._matches) * 2
    residuals = flex.double(nelem)
    jacobian_t = flex.double(flex.grid(
        len(self._prediction_parameterisation), nelem))
    weights = flex.double(nelem)

    for i, m in enumerate(self._matches):
      residuals[2*i] = m.x_resid
      residuals[2*i + 1] = m.y_resid
      #residuals[3*i + 2] = m.Phiresid

      # are these the right weights? Or inverse, or sqrt?
      weights[2*i] = m.weight_x_obs
      weights[2*i + 1] = m.weight_y_obs
      #weights[3*i + 2] = m.weightPhio

      # m.gradients is a nparam length list, each element of which is a
      # doublet of values, (dX/dp_n, dY/dp_n)
      dX_dp, dY_dp = zip(*m.gradients)

      # FIXME Here we paste columns into the Jacobian transpose then take
      # its transpose when done. This seems inefficient: can we just start
      # with the Jacobian and fill elements sequentially, using row-major
      # order to ensure the values are filled in the right order?

      # fill jacobian elements here.
      jacobian_t.matrix_paste_column_in_place(flex.double(dX_dp), 2*i)
      jacobian_t.matrix_paste_column_in_place(flex.double(dY_dp), 2*i+1)
      #jacobian_t.matrix_paste_column_in_place(flex.double(dPhi_dp), 3*i+2)

    # We return the Jacobian, not its transpose.
    jacobian_t.matrix_transpose_in_place()

    return(residuals, jacobian_t, weights)

  def compute_functional_and_gradients(self):
    """calculate the value of the target function and its gradients"""

    self._matches = self._reflection_manager.get_matches()
    self._nref = len(self._matches)

    # This is a hack for the case where nref=0. This should not be necessary
    # if bounds are provided for parameters to stop the algorithm exploring
    # unreasonable regions of parameter space where no predictions exist.
    # Unfortunately the L-BFGS line search does make such extreme trials.
    if self._nref == 0:
      return 1.e12, [1.] * len(self._prediction_parameterisation)

    # compute target function
    L = 0.5 * sum([m.weight_x_obs * m.x_resid2 +
                   m.weight_y_obs * m.y_resid2
                   for m in self._matches])

    # prepare list of gradients
    dL_dp = [0.] * len(self._prediction_parameterisation)

    # the gradients wrt each parameter are stored with the matches
    for m in self._matches:

      for j, (grad_X, grad_Y) in enumerate(m.gradients):
        dL_dp[j] += (m.weight_x_obs * m.x_resid * grad_X +
                     m.weight_y_obs * m.y_resid * grad_Y)

    return (L, dL_dp)

  def curvatures(self):
    """First order approximation to the diagonal of the Hessian based on the
    least squares form of the target"""

    # This is a hack for the case where nref=0. This should not be necessary
    # if bounds are provided for parameters to stop the algorithm exploring
    # unreasonable regions of parameter space where there are no predictions
    if self._nref == 0:
      return [1.] * len(self._prediction_parameterisation)

    # prepare lists of gradients and curvatures
    curv = [0.] * len(self._prediction_parameterisation)

    # for each reflection, get the approximate curvatures wrt each parameter
    for m in self._matches:

      for j, (grad_X, grad_Y) in enumerate(m.gradients):
        curv[j] += (m.weight_x_obs * grad_X**2 +
                    m.weight_y_obs * grad_Y**2)

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert all([c > 0.0 for c in curv])

    return curv

  def rmsds(self):
    """calculate unweighted RMSDs"""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()

    resid_x = sum((m.x_resid2 for m in self._matches))
    resid_y = sum((m.y_resid2 for m in self._matches))

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

class ReflectionManagerXY(ReflectionManager):
  """Overloads for a Reflection Manager that does not exclude
  reflections too close to the spindle, and reports only information
  about X, Y residuals"""

  def _spindle_beam_plane_normal(self):
    """There is no goniometer, so overload to return None"""

    return None

  def _id_refs_to_keep(self, obs_data):
    """For this version of the class, only reject the (0,0,0) reflections.
    We don't want to exclude reflections close to the spindle, as the spindle
    may not exist"""

    inc = [i for i, ref in enumerate(obs_data) if ref.miller_index != (0,0,0)]

    return inc

  def _create_working_set(self, indices):
    """Make a subset of the indices of reflections to use in refinement.

    This version ignores nref_per_degree"""

    working_indices = indices
    sample_size = len(working_indices)

    # set maximum sample size
    if self._max_num_obs:
      if sample_size > self._max_num_obs:
        sample_size = self._max_num_obs

    # sample the data and record the sample size
    if sample_size < len(working_indices):
      self._sample_size = sample_size
      working_indices = random.sample(working_indices,
                                      self._sample_size)
    return(working_indices)

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()

    if self._verbosity > 1:

      from scitbx.math import five_number_summary
      x_resid = [e.x_resid for e in l]
      y_resid = [e.y_resid for e in l]
      w_x = [e.weight_x_obs for e in l]
      w_y = [e.weight_y_obs for e in l]

      print "\nSummary statistics for observations matched to predictions:"
      print ("                      "
             "Min         Q1        Med         Q3        Max")
      print "(Xc-Xo)        {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(x_resid))
      print "(Yc-Yo)        {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(y_resid))
      print "X weights      {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_x))
      print "Y weights      {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_y))
      print

      if len(l) >= 20 and self._verbosity > 2:

        sl = self._sort_obs_by_residual(l)
        print "Reflections with the worst 20 positional residuals:"
        print "H, K, L, x_resid, y_resid, weight_x_obs, weight_y_obs"
        fmt = "(%3d, %3d, %3d) %5.3f %5.3f %5.3f %5.3f"
        for i in xrange(20):
          e = sl[i]
          msg = fmt % tuple(e.miller_index + (e.x_resid,
                           e.y_resid,
                           e.weight_x_obs,
                           e.weight_y_obs))
          print msg
        print

    return

  def reject_outliers(self):
    """Unset the use flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    matches = [obs for obs in self._obs_pred_pairs if obs.is_matched]

    x_resid = [e.x_resid for e in matches]
    y_resid = [e.y_resid for e in matches]

    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
    min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)

    iqr_x = q3_x - q1_x
    iqr_y = q3_y - q1_y

    cut_x = self._iqr_multiplier * iqr_x
    cut_y = self._iqr_multiplier * iqr_y

    for m in matches:
      if m.x_resid > q3_x + cut_x: m.is_matched = False
      if m.x_resid < q1_x - cut_x: m.is_matched = False
      if m.y_resid > q3_y + cut_y: m.is_matched = False
      if m.y_resid < q1_y - cut_y: m.is_matched = False

    nreject = [m.is_matched for m in matches].count(False)

    if nreject == 0: return False

    if self._verbosity > 1:
      print "%d reflections have been rejected as outliers" % nreject

    return True
