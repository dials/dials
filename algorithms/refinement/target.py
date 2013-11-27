#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
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
from math import pi, sqrt
from cctbx.array_family import flex
import random

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

  rmsd_names = ["RMSD_X", "RMSD_Y", "RMSD_Phi"]

  def __init__(self, reflection_predictor, detector,
               ref_manager,
               prediction_parameterisation):

    self._reflection_predictor = reflection_predictor
    self._detector = detector
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
    self._prediction_parameterisation.prepare()

    # reset the 'use' flag for all observations
    self._reflection_manager.reset_accepted_reflections()

    # loop over all reflections in the manager
    for obs in self._reflection_manager.get_obs():

      # get data from the observation
      h = obs.miller_index
      frame = obs.frame_obs
      panel = obs.panel

      # duck-typing for scan varying version of
      # prediction_parameterisation
      try:

        # compose the prediction parameterisation at the
        # requested image number
        self._prediction_parameterisation.compose(frame)

        # extract UB matrix
        UB = self._prediction_parameterisation.get_UB(frame)

        # predict for this hkl at setting UB
        predictions = self._reflection_predictor.predict(h, UB)

      except AttributeError:

        # predict for this hkl
        predictions = self._reflection_predictor.predict(h)

      # obtain the impact positions
      impacts = ray_intersection(self._detector, predictions,
                                 panel=panel)

      # find the prediction with the right 'entering' flag
      try:
        i = [x.entering == obs.entering \
             for x in impacts].index(True)
      except ValueError:
        # we don't have a prediction for this obs
        continue

      ref = impacts[i]
      x_calc, y_calc = ref.image_coord_mm

      # do not wrap around multiples of 2*pi; keep the full rotation
      # from zero to differentiate repeat observations.
      resid = ref.rotation_angle - (obs.phi_obs % TWO_PI)

      # ensure this is the smaller of two possibilities
      resid = (resid + pi) % TWO_PI - pi

      phi_calc = obs.phi_obs + resid
      s_calc = matrix.col(ref.beam_vector)

      # calculate gradients for this reflection
      grads = self._prediction_parameterisation.get_gradients(
                                  h, s_calc, phi_calc, panel, frame)

      # store all this information in the matched obs-pred pair
      obs.update_prediction(x_calc, y_calc, phi_calc, s_calc, grads)

    if self._reflection_manager.first_update:

      # delete all obs-pred pairs from the manager that do not
      # have a prediction
      self._reflection_manager.strip_unmatched_observations()

      self._reflection_manager.first_update = False

    return

  def get_num_matches(self):
    """return the number of reflections currently used in the calculation"""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()

    return len(self._matches)

  def compute_functional_and_gradients(self):
    """calculate the target function value and its gradients"""

    # To be implemented by a derived class
    raise RuntimeError('implement me')

  def compute_residuals_and_gradients(self):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    # To be implemented by a derived class
    raise RuntimeError('implement me')

  def curvatures(self):
    """First order approximation to the diagonal of the Hessian based on the
    least squares form of the target"""

    # To be implemented by a derived class
    raise RuntimeError('implement me')

  def rmsds(self):
    """calculate unweighted RMSDs"""

    # To be implemented by a derived class
    raise RuntimeError('implement me')

  def achieved(self):
    """return True to terminate the refinement. To be implemented by
    a derived class"""

    # To be implemented by a derived class
    raise RuntimeError('implement me')

class LeastSquaresPositionalResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y and phi, terminating on achieved
  rmsd (or on intrisic convergence of the chosen minimiser)"""

  def __init__(self, reflection_predictor, detector, ref_man,
               prediction_parameterisation,
               image_width, frac_binsize_cutoff=0.33333,
               absolute_cutoffs=None):

    Target.__init__(self, reflection_predictor, detector, ref_man,
                    prediction_parameterisation)

    # Set up the RMSD achieved criterion
    if not absolute_cutoffs:
      pixel_sizes = [p.get_pixel_size() for p in detector]
      min_px_size_x = min(e[0] for e in pixel_sizes)
      min_px_size_y = min(e[1] for e in pixel_sizes)
      self._binsize_cutoffs = [min_px_size_x * frac_binsize_cutoff,
                               min_px_size_y * frac_binsize_cutoff,
                               image_width * frac_binsize_cutoff]
    else:
      assert len(absolute_cutoffs) == 3
      self._binsize_cutoffs = absolute_cutoffs

    return

  def compute_residuals_and_gradients(self):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    self._matches = self._reflection_manager.get_matches()

    # return residuals and weights as 1d flex.double vectors
    nelem = len(self._matches) * 3
    residuals = flex.double(nelem)
    jacobian_t = flex.double(flex.grid(
        len(self._prediction_parameterisation), nelem))
    weights = flex.double(nelem)

    for i, m in enumerate(self._matches):
      residuals[3*i] = m.x_resid
      residuals[3*i + 1] = m.y_resid
      residuals[3*i + 2] = m.phi_resid

      # are these the right weights? Or inverse, or sqrt?
      weights[3*i] = m.weight_x_obs
      weights[3*i + 1] = m.weight_y_obs
      weights[3*i + 2] = m.weight_phi_obs

      # m.gradients is a nparam length list, each element of which is a
      # triplet of values, (dX/dp_n, dY/dp_n, dPhi/dp_n)
      dX_dp, dY_dp, dPhi_dp = zip(*m.gradients)

      # FIXME Here we paste columns into the Jacobian transpose then take
      # its transpose when done. This seems inefficient: can we just start
      # with the Jacobian and fill elements sequentially, using row-major
      # order to ensure the values are filled in the right order?

      # fill jacobian elements here.
      jacobian_t.matrix_paste_column_in_place(flex.double(dX_dp), 3*i)
      jacobian_t.matrix_paste_column_in_place(flex.double(dY_dp), 3*i+1)
      jacobian_t.matrix_paste_column_in_place(flex.double(dPhi_dp), 3*i+2)

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

      for j, (grad_X, grad_Y, grad_Phi) in enumerate(m.gradients):
        curv[j] += (m.weight_x_obs * grad_X**2 +
                    m.weight_y_obs * grad_Y**2 +
                    m.weight_phi_obs * grad_Phi**2)

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert all([c > 0.0 for c in curv])

    return curv

  def rmsds(self):
    """calculate unweighted RMSDs"""

    if not self._matches:
      self._matches = self._reflection_manager.get_matches()

    resid_x = sum((m.x_resid2 for m in self._matches))
    resid_y = sum((m.y_resid2 for m in self._matches))
    resid_phi = sum((m.phi_resid2 for m in self._matches))

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

class ObsPredMatch:
  """
  A bucket class containing data for a prediction that has been
  matched to an observation.

  This contains all the raw material needed to calculate the target function
  value and gradients
  """

  # initialise with an observation
  def __init__(self, iobs, crystal, hkl, entering, frame, panel,
                     x_obs, sigx_obs, weight_x_obs,
                     y_obs, sigy_obs, weight_y_obs,
                     phi_obs, sigphi_obs, weight_phi_obs):

    self.iobs = iobs
    self.crystal_id = crystal
    self.miller_index = hkl
    self.entering = entering
    self.frame_obs = frame
    self.panel = panel
    self.x_obs = x_obs
    self.sigx_obs = sigx_obs
    self.weight_x_obs = weight_x_obs

    self.y_obs = y_obs
    self.sigy_obs = sigy_obs
    self.weight_y_obs = weight_y_obs

    self.phi_obs = phi_obs
    self.sigphi_obs = sigphi_obs
    self.weight_phi_obs = weight_phi_obs

    self.x_calc = None
    self.y_calc = None
    self.phi_calc = None
    self.s_calc = None

    # gradients will be a list, of length equal to the number of free
    # parameters, whose elements are the gradients in the space of
    # the residuals, e.g. (dX/dp, dY/dp, dPhi/dp)
    self.gradients = None

    self.y_resid = None
    self.y_resid2 = None
    self.x_resid = None
    self.x_resid2 = None
    self.phi_resid = None
    self.phi_resid2 = None

    self.is_matched = False

  # update with a prediction
  def update_prediction(self, x_calc, y_calc, phi_calc, s_calc, gradients):

    self.x_calc = x_calc
    self.y_calc = y_calc
    self.phi_calc = phi_calc
    self.s_calc = s_calc

    self.gradients = gradients

    # calculate residuals
    self.x_resid = x_calc - self.x_obs
    self.x_resid2 = self.x_resid**2
    self.y_resid = y_calc - self.y_obs
    self.y_resid2 = self.y_resid**2
    self.phi_resid = phi_calc - self.phi_obs
    self.phi_resid2 = self.phi_resid**2

    self.is_matched = True

  def reset(self):

    """Flag this observation to not be used"""
    self.is_matched = False


class ReflectionManager(object):
  """A class to maintain information about observed and predicted
  reflections for refinement."""

  def __init__(self, reflections,
                     beam, gonio,
                     sweep_range_rad=None,
                     nref_per_degree=None,
                     min_num_obs=20,
                     max_num_obs=None,
                     inclusion_cutoff=0.1,
                     verbosity=0):

    # track whether this is the first update of predictions or not
    self.first_update = True

    # set verbosity
    self._verbosity = verbosity

    # keep references to the beam, goniometer and sweep range (for
    # reflection exclusion and subsetting)
    self._beam = beam
    self._gonio = gonio
    self._sweep_range_rad = sweep_range_rad

    # find vector normal to the spindle-beam plane for the initial model
    self._vecn = self._spindle_beam_plane_normal()

    # set up the reflection inclusion cutoff
    self._inclusion_cutoff = inclusion_cutoff

    # exclude reflections that fail inclusion criteria
    self._input_size = len(reflections)
    refs_to_keep = self._id_refs_to_keep(reflections)
    self._accepted_refs_size = len(refs_to_keep)

    # choose a random subset of data for refinement
    self._sample_size = self._accepted_refs_size
    self._nref_per_degree = nref_per_degree
    self._max_num_obs = max_num_obs
    refs_to_keep = self._create_working_set(refs_to_keep)

    # store observation information in a list of observation-prediction
    # pairs (prediction information will go in here later)
    self._obs_pred_pairs = []
    for i in refs_to_keep:

      ref = reflections[i]
      crystal = ref.crystal
      h = ref.miller_index
      s = matrix.col(ref.beam_vector)
      entering = s.dot(self._vecn) < 0.
      frame = ref.frame_number
      panel = ref.panel_number
      x = ref.centroid_position[0]
      y = ref.centroid_position[1]
      phi = ref.centroid_position[2]
      sig_x, sig_y, sig_phi = [sqrt(e) for e in ref.centroid_variance]
      w_x, w_y, w_phi = [1. / e for e in ref.centroid_variance]

      self._obs_pred_pairs.append(ObsPredMatch(i, crystal, h,
                                               entering, frame, panel,
                                               x, sig_x, w_x,
                                               y, sig_y, w_y,
                                               phi, sig_phi, w_phi))

    # fail if there are too few reflections in the manager
    self._min_num_obs = min_num_obs
    if len(self._obs_pred_pairs) < self._min_num_obs:
      msg = ('Remaining number of reflections = {0}, which is below '+ \
          'the configured limit for creating this reflection ' + \
          'manager').format(len(self._obs_pred_pairs))
      raise RuntimeError(msg)

  def _spindle_beam_plane_normal(self):
    """return a unit vector that when placed at the origin of reciprocal
    space, points to the hemisphere of the Ewald sphere
    in which reflections cross from inside to outside of the sphere"""

    # NB vector in +ve Y direction when using imgCIF coordinate frame
    return matrix.col(self._beam.get_s0()).cross(
                    matrix.col(self._gonio.get_rotation_axis())).normalize()

  def _id_refs_to_keep(self, obs_data):
    """Create a selection of observations that pass certain conditions.

    This includes outlier rejection plus rejection of reflections
    too close to the spindle, and rejection of the (0,0,0) Miller index"""

    # TODO Add outlier rejection. Should use robust statistics.
    # See notes on M-estimators from Garib (when I have them)

    axis = matrix.col(self._gonio.get_rotation_axis())
    s0 = matrix.col(self._beam.get_s0())

    inc = [i for i, ref in enumerate(obs_data) if ref.miller_index != (0,0,0) \
           and self._inclusion_test(matrix.col(ref.beam_vector), axis, s0)]

    return inc

  def _inclusion_test(self, s, axis, s0):
    """Test scattering vector s for inclusion"""

    # reject reflections for which the parallelepiped formed between
    # the gonio axis, s0 and s has a volume of less than the cutoff.
    # Those reflections are by definition closer to the spindle-beam
    # plane and for low values of the cutoff are troublesome to
    # integrate anyway.

    test = abs(axis.dot(matrix.col(s).cross(s0))) > \
        self._inclusion_cutoff

    return test

  def _create_working_set(self, indices):
    """Make a subset of the indices of reflections to use in refinement"""

    working_indices = indices
    sample_size = len(working_indices)

    # set sample size according to nref_per_degree
    if self._nref_per_degree:
      width = abs(self._sweep_range_rad[1] -
                  self._sweep_range_rad[0]) * RAD_TO_DEG
      sample_size = int(self._nref_per_degree * width)

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

  def get_input_size(self):
    """Return the number of observations in the initial list supplied
    as input"""

    return self._input_size

  def get_accepted_refs_size(self):
    """Return the number of observations that pass inclusion criteria and
    can potentially be used for refinement"""

    return self._accepted_refs_size

  def get_sample_size(self):
    """Return the number of observations in the working set to be
    used for refinement"""

    return self._sample_size

  def _sort_obs_by_residual(self, obs, angular=False):
    """For diagnostic purposes, sort the obs-pred matches so that the
    highest residuals are first. By default, sort by positional
    residual, unless angular=True.

    The earliest entries in the return list may be those that are
    causing problems in refinement.

    """

    if angular:
      k = lambda e: e.phi_resid
    else:
      k = lambda e: e.x_resid**2 + e.y_resid**2
    sort_obs = sorted(obs, key=k, reverse=True)

    return sort_obs

  def get_matches(self, silent = False):
    """For every observation matched with a prediction return all data"""

    l = [obs for obs in self._obs_pred_pairs if obs.is_matched]

    if self._verbosity > 2 and len(l) > 20 and not silent:

      sl = self._sort_obs_by_residual(l)
      print "Reflections with the worst 20 positional residuals:"
      print "H, K, L, x_resid, y_resid, phi_resid, weight_x_obs, weight_y_obs, " + \
            "weight_phi_obs"
      fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
      for i in xrange(20):
        e = sl[i]
        msg = fmt % tuple(e.miller_index + (e.x_resid,
                         e.y_resid,
                         e.phi_resid,
                         e.weight_x_obs,
                         e.weight_y_obs,
                         e.weight_phi_obs))
        print msg
      print
      sl = self._sort_obs_by_residual(l, angular=True)
      print "\nReflections with the worst 20 angular residuals:"
      print "H, K, L, x_resid, y_resid, phi_resid, weight_x_obs, weight_y_obs, " + \
            "weight_phi_obs"
      fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
      for i in xrange(20):
        e = sl[i]
        msg = fmt % tuple(e.miller_index + (e.x_resid,
                         e.y_resid,
                         e.phi_resid,
                         e.weight_x_obs,
                         e.weight_y_obs,
                         e.weight_phi_obs))
        print msg
      print

    return l

  def strip_unmatched_observations(self):
    """Delete observations from the manager that are not matched to a
    prediction. Typically used once, after the first update of
    predictions."""

    self._obs_pred_pairs = [e for e in self._obs_pred_pairs if e.is_matched]

    if len(self._obs_pred_pairs) < self._min_num_obs:
      msg = ('Remaining number of reflections = {0}, which is below '+ \
          'the configured limit for this reflection manager').format(
              len(self._obs_pred_pairs))
      raise RuntimeError(msg)

    if self._verbosity > 1:
      print len(self._obs_pred_pairs), "reflections remain in the manager after " + \
          "removing those unmatched with predictions"

    return

  def reset_accepted_reflections(self):
    """Reset all observations to use=False in preparation for a new set of
    predictions"""

    for obs in self._obs_pred_pairs: obs.reset()

  def get_obs(self):
    """Get the list of managed observations"""

    return self._obs_pred_pairs
