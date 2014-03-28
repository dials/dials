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

# dials imports
from dials.algorithms.spot_prediction import ray_intersection

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
               prediction_parameterisation):

    self._reflection_predictor = reflection_predictor
    self._experiments = experiments
    self._reflection_manager = ref_manager
    self._prediction_parameterisation = prediction_parameterisation

    # Quantities to cache each step
    self._rmsds = None
    self._matches = None

    return

  def predict(self):
    """perform reflection prediction and update the reflection manager"""

    # update the reflection_predictor and the prediction parameterisation
    # with the scan-independent part of the current geometry
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
    resid = phi_calc - (flex.fmod(phi_obs, TWO_PI))
    # ensure this is the smaller of two possibilities
    resid = flex.fmod((resid + pi), TWO_PI) - pi
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
    self._matches = self._reflection_manager.get_matches()

    return

  def calculate_gradients(self):
    """delegate to the prediction_parameterisation object to calculate
    gradients for all the matched reflections."""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()
    self._gradients = self._prediction_parameterisation.get_gradients(
      self._matches)

    return self._gradients

  def get_num_matches(self):
    """return the number of reflections currently used in the calculation"""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()

    return len(self._matches)

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
               absolute_cutoffs=None):

    Target.__init__(self, experiments, reflection_predictor, ref_man,
                    prediction_parameterisation)

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

  def compute_residuals_and_gradients(self):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    dX_dp, dY_dp, dPhi_dp = self.calculate_gradients()

    # return residuals and weights as 1d flex.double vectors
    nelem = len(self._matches) * 3
    nparam = len(self._prediction_parameterisation)
    residuals = flex.double.concatenate(self._matches['x_resid'],
                                        self._matches['y_resid'])
    residuals.extend(self._matches['phi_resid'])
    #jacobian_t = flex.double(flex.grid(
    #    len(self._prediction_parameterisation), nelem))
    weights, w_y, w_z = self._matches['xyzobs.mm.weights'].parts()
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

  # FIXME
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
                   m.weight_y_obs * m.y_resid2 +
                   m.weight_phi_obs * m.phi_resid2
                   for m in self._matches])

    # prepare list of gradients
    dL_dp = [0.] * len(self._prediction_parameterisation)

    # the gradients wrt each parameter are stored with the matches
    for m in self._matches:

      for j, (grad_X, grad_Y, grad_Phi) in enumerate(m.gradients):
        dL_dp[j] += (m.weight_x_obs * m.x_resid * grad_X +
                     m.weight_y_obs * m.y_resid * grad_Y +
                     m.weight_phi_obs * m.phi_resid * grad_Phi)

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

      for j, (grad_x, grad_y, grad_phi) in enumerate(m.gradients):
        curv[j] += (m.weight_x_obs * grad_x**2 +
                    m.weight_y_obs * grad_y**2 +
                    m.weight_phi_obs * grad_phi**2)

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert all([c > 0.0 for c in curv])

    return curv

  def rmsds(self):
    """calculate unweighted RMSDs"""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()

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

class LeastSquaresXYResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y only, terminating on achieved
  rmsd (or on intrisic convergence of the chosen minimiser)"""

  rmsd_names = ["RMSD_X", "RMSD_Y"]

  def __init__(self, experiments, reflection_predictor, ref_man,
               prediction_parameterisation,
               frac_binsize_cutoff=0.33333,
               absolute_cutoffs=None):

    Target.__init__(self, experiments, reflection_predictor, ref_man,
                    prediction_parameterisation)

    # Set up the RMSD achieved criterion. For simplicity, we take detector from
    # the first Experiment only.
    detector = experiments[0].detector
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

