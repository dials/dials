#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Contains classes used to construct a target function for refinement,
principally Target and ReflectionManager."""

# python and cctbx imports
from __future__ import division
from scitbx import matrix
from math import pi, sqrt, ceil
from cctbx.array_family import flex
import random
import abc

# constants
TWO_PI = 2.0 * pi
RAD_TO_DEG = 180. / pi

class Target(object):
  """Abstract interface for a target function class

  A Target object will be used by a Refinery. It will refer to a Reflection
  Manager to get a list of observations. It will perform reflection prediction
  on those observations and update the reflection manager with those
  predictions. It will then query the reflection manager to get a list of
  accepted reflections with observed and calculated positions. These are the
  reflections for use in refinement. It obtains the gradients of reflection
  positions from a relevant prediction parameterisation object. With all of
  this information in place, it calculates the value of the target function,
  the gradients of the target function and auxiliary information (e.g. RMSDs).

  Concrete instances of this class implement the actual target function
  calculation. The base class should not determine what type of target
  function is used (e.g. least squares target), or limit whether the object
  is used for a detector space & phi residual, or a reciprocal space residual.
  This should all be set by a derived class.
  """

  __metaclass__  = abc.ABCMeta
  rmsd_names = ["RMSD_X", "RMSD_Y", "RMSD_Phi"]

  def __init__(self, experiments, reflection_predictor, ref_manager,
               prediction_parameterisation, jacobian_max_nref=None):

    self._reflection_predictor = reflection_predictor
    self._experiments = experiments
    self._reflection_manager = ref_manager
    self._prediction_parameterisation = prediction_parameterisation

    # Quantities to cache each step
    self._rmsds = None
    self._matches = None

    # Keep maximum number of reflections used to for Jacobian calculation, if
    # a cutoff is required
    self._jacobian_max_nref = jacobian_max_nref
    self._finished_residuals_and_gradients = False

    return

  def predict(self):
    """perform reflection prediction and update the reflection manager"""

    # update the reflection_predictor with the scan-independent part of the
    # current geometry
    self._reflection_predictor.update()

    # reset the 'use' flag for all observations
    self._reflection_manager.reset_accepted_reflections()

    # duck-typing for VaryingCrystalPredictionParameterisation. Only this
    # class has a compose(reflections) method. Sets ub_matrix (and caches
    # derivatives).
    try:
      # TODO set U, B and UB arrays in the reflection_table. Cache derivatives of these
      # in the _prediction_parameterisation object (so we don't have to compose
      # at each scan point again later for the derivatives). If this is the
      # first prediction of reflections some reflections could then be rejected
      # by the reflection manager, and
      # then it would have been a waste of time calculating derivatives for
      # them. However, it is expected that scan-varying refinement is done after
      # scan static refinement, so probably all of the rejections have been done
      # already (so don't worry, be happy)
      self._prediction_parameterisation.compose(
          self._reflection_manager.get_obs())
    except AttributeError:
      pass

    # do prediction (updates reflection table in situ). Scan-varying prediction
    # is done automatically if the crystal has scan-points (assuming reflections
    # have ub_matrix set)
    reflections = self._reflection_manager.get_obs()
    self._reflection_predictor.predict(reflections)

    x_obs, y_obs, phi_obs = reflections['xyzobs.mm.value'].parts()
    x_calc, y_calc, phi_calc = reflections['xyzcal.mm'].parts()
    # do not wrap around multiples of 2*pi; keep the full rotation
    # from zero to differentiate repeat observations.
    resid = phi_calc - (flex.fmod_positive(phi_obs, TWO_PI))
    # ensure this is the smaller of two possibilities
    resid = flex.fmod_positive((resid + pi), TWO_PI) - pi
    phi_calc = phi_obs + resid
    # put back in the reflections
    reflections['xyzcal.mm'] = flex.vec3_double(x_calc, y_calc, phi_calc)

    # calculate residuals and assign columns
    reflections['x_resid'] = x_calc - x_obs
    reflections['x_resid2'] = reflections['x_resid']**2
    reflections['y_resid'] = y_calc - y_obs
    reflections['y_resid2'] = reflections['y_resid']**2
    reflections['phi_resid'] = phi_calc - phi_obs
    reflections['phi_resid2'] = reflections['phi_resid']**2

    # set used_in_refinement flag to all those that had predictions
    mask = reflections.get_flags(reflections.flags.predicted)
    reflections.set_flags(mask, reflections.flags.used_in_refinement)

    # collect the matches
    self.update_matches(force=True)

    return

  def predict_for_reflection_table(self, reflections):
    """perform prediction for all reflections in the supplied table"""

    # ensure the predictor is ready with the right models
    self._reflection_predictor.update()
    # set the entering flags
    from dials.algorithms.refinement.reflection_manager import calculate_entering_flags
    reflections['entering'] = calculate_entering_flags(reflections, self._experiments)
    self._reflection_predictor.predict(reflections)
    x_obs, y_obs, phi_obs = reflections['xyzobs.mm.value'].parts()
    x_calc, y_calc, phi_calc = reflections['xyzcal.mm'].parts()
    # do not wrap around multiples of 2*pi; keep the full rotation
    # from zero to differentiate repeat observations.
    resid = phi_calc - (flex.fmod(phi_obs, TWO_PI))
    # ensure this is the smaller of two possibilities
    resid = flex.fmod((resid + pi), TWO_PI) - pi
    phi_calc = phi_obs + resid
    # put back in the reflections
    reflections['xyzcal.mm'] = flex.vec3_double(x_calc, y_calc, phi_calc)
    return reflections

  def calculate_gradients(self, reflections=None):
    """delegate to the prediction_parameterisation object to calculate
    gradients for all the matched reflections, or just for those specified"""

    self.update_matches()

    if reflections:
      self._gradients = self._prediction_parameterisation.get_gradients(
        reflections)
    else:
      self._gradients = self._prediction_parameterisation.get_gradients(
        self._matches)

    return self._gradients

  def get_num_matches(self):
    """return the number of reflections currently used in the calculation"""

    self.update_matches()
    return len(self._matches)

  def update_matches(self, force=False):
    """ensure the observations matched to predictions are up to date"""

    if not self._matches or force:
      self._matches = self._reflection_manager.get_matches()

    return

  @property
  def finished_residuals_and_gradients(self):
    """used when looping through subsets of the full jacobian to determine when
    to break the loop. Reset flag whenever called"""

    val = self._finished_residuals_and_gradients
    self._finished_residuals_and_gradients = False
    return val

  @abc.abstractmethod
  def compute_functional_and_gradients(self):
    """calculate the target function value and its gradients"""

    pass

  @abc.abstractmethod
  def compute_residuals_and_gradients(self):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    pass

  @abc.abstractmethod
  def curvatures(self):
    """First order approximation to the diagonal of the Hessian based on the
    least squares form of the target"""

    pass

  @abc.abstractmethod
  def rmsds(self):
    """calculate unweighted RMSDs"""

    pass

  @abc.abstractmethod
  def achieved(self):
    """return True to terminate the refinement."""

    pass

class LeastSquaresPositionalResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y and phi, terminating on achieved
  rmsd (or on intrisic convergence of the chosen minimiser)"""

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
    if experiments[0].scan:
      temp = experiments[0].scan.get_oscillation(deg=False)
      image_width_rad = temp[1] - temp[0]
    else:
      image_width_rad = None
    detector = experiments[0].detector

    if not absolute_cutoffs:
      pixel_sizes = [p.get_pixel_size() for p in detector]
      min_px_size_x = min(e[0] for e in pixel_sizes)
      min_px_size_y = min(e[1] for e in pixel_sizes)
      self._binsize_cutoffs = [min_px_size_x * frac_binsize_cutoff,
                               min_px_size_y * frac_binsize_cutoff,
                               image_width_rad * frac_binsize_cutoff]
    else:
      assert len(absolute_cutoffs) == 3
      self._binsize_cutoffs = absolute_cutoffs

    # predict reflections and finalise reflection manager
    self.predict()
    self._reflection_manager.finalise()

    return

  def compute_residuals_and_gradients(self, block_num=0):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    self.update_matches()
    if self._jacobian_max_nref:
      start = block_num * self._jacobian_max_nref
      end = (block_num + 1) * self._jacobian_max_nref
      matches = self._matches[start:end]
      self._finished_residuals_and_gradients = True if \
        end >= len(self._matches) else False
    else:
      matches = self._matches
      self._finished_residuals_and_gradients = True

    dX_dp, dY_dp, dPhi_dp = self.calculate_gradients(matches)

    # return residuals and weights as 1d flex.double vectors
    nelem = len(matches) * 3
    nparam = len(self._prediction_parameterisation)
    residuals = flex.double.concatenate(matches['x_resid'],
                                        matches['y_resid'])
    residuals.extend(matches['phi_resid'])
    #jacobian_t = flex.double(flex.grid(
    #    len(self._prediction_parameterisation), nelem))
    weights, w_y, w_z = matches['xyzobs.mm.weights'].parts()
    weights.extend(w_y)
    weights.extend(w_z)

    jacobian = flex.double(flex.grid(nelem, nparam))
    # loop over parameters
    for i in range(nparam):
      dX, dY, dPhi = dX_dp[i], dY_dp[i], dPhi_dp[i]
      col = flex.double.concatenate(dX, dY)
      col.extend(dPhi)
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
    phi_resid = self._matches['phi_resid']
    x_resid2 = self._matches['x_resid2']
    y_resid2 = self._matches['y_resid2']
    phi_resid2 = self._matches['phi_resid2']
    w_x, w_y, w_phi = self._matches['xyzobs.mm.weights'].parts()

    # calculate target function
    temp = w_x * x_resid2 + w_y * y_resid2 + w_phi * phi_resid2
    L = 0.5 * flex.sum(temp)

    # prepare list of gradients
    dL_dp = [0.] * len(self._prediction_parameterisation)

    dX_dp, dY_dp, dPhi_dp = self.calculate_gradients()

    w_x_resid = w_x * x_resid
    w_y_resid = w_y * y_resid
    w_phi_resid = w_phi * phi_resid

    for i in range(len(self._prediction_parameterisation)):
      dX, dY, dPhi = dX_dp[i], dY_dp[i], dPhi_dp[i]
      temp = w_x_resid * dX + w_y_resid * dY + w_phi_resid * dPhi
      dL_dp[i] = flex.sum(temp)

    return (L, dL_dp)

  def curvatures(self):
    """First order approximation to the diagonal of the Hessian based on the
    least squares form of the target"""

    # relies on compute_functional_and_gradients being called first
    dX_dp, dY_dp, dPhi_dp = self._gradients
    w_x, w_y, w_phi = self._matches['xyzobs.mm.weights'].parts()

    # This is a hack for the case where nref=0. This should not be necessary
    # if bounds are provided for parameters to stop the algorithm exploring
    # unreasonable regions of parameter space where there are no predictions
    if self._nref == 0:
      return [1.] * len(self._prediction_parameterisation)

    # prepare lists of gradients and curvatures
    curv = [0.] * len(self._prediction_parameterisation)

    # for each reflection, get the approximate curvatures wrt each parameter
    for i in range(len(self._prediction_parameterisation)):
      dX, dY, dPhi = dX_dp[i], dY_dp[i], dPhi_dp[i]
      temp = w_x * dX**2 + w_y * dY**2 + w_phi * dPhi**2
      curv[i] = flex.sum(temp)

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert all([c > 0.0 for c in curv])

    return curv

  def rmsds(self):
    """calculate unweighted RMSDs"""

    self.update_matches()
    resid_x = flex.sum(self._matches['x_resid2'])
    resid_y = flex.sum(self._matches['y_resid2'])
    resid_phi = flex.sum(self._matches['phi_resid2'])

    # cache rmsd calculation for achieved test
    n = len(self._matches)
    self._rmsds = (sqrt(resid_x / n),
                   sqrt(resid_y / n),
                   sqrt(resid_phi / n))

    return self._rmsds

  def achieved(self):
    """RMSD criterion for target achieved """
    r = self._rmsds if self._rmsds else self.rmsds()

    # reset cached rmsds to avoid getting out of step
    self._rmsds = None

    if (r[0] < self._binsize_cutoffs[0] and
        r[1] < self._binsize_cutoffs[1] and
        r[2] < self._binsize_cutoffs[2]):
      return True
    return False
