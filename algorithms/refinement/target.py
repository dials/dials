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
from math import pi, sqrt, ceil
from cctbx.array_family import flex
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
  _grad_names = ['dX_dp', 'dY_dp', 'dphi_dp']
  rmsd_names = ["RMSD_X", "RMSD_Y", "RMSD_Phi"]
  rmsd_units = ["mm", "mm", "rad"]

  def __init__(self, experiments, reflection_predictor, ref_manager,
               prediction_parameterisation, restraints_parameterisation=None,
               gradient_calculation_blocksize=None):

    self._reflection_predictor = reflection_predictor
    self._experiments = experiments
    self._reflection_manager = ref_manager
    self._prediction_parameterisation = prediction_parameterisation
    self._restraints_parameterisation = restraints_parameterisation

    # Quantities to cache each step
    self._rmsds = None
    self._matches = None

    # Keep maximum number of reflections used for Jacobian calculation, if
    # a cutoff is required
    self._gradient_calculation_blocksize = gradient_calculation_blocksize

    return

  def set_prediction_parameterisation(self, prediction_parameterisation):
    """For circumstances where the PredictionParameterisation object was not
    available at initialisation, set it with this method"""

    self._prediction_parameterisation = prediction_parameterisation
    return

  def set_restraints_parameterisation(self, restraints_parameterisation):
    """For circumstances where the RestraintsParameterisation object was not
    available at initialisation, set it with this method"""

    self._restraints_parameterisation = restraints_parameterisation
    return

  def _predict_core(self, reflections):
    """perform prediction for the specified reflections"""

    # update the reflection_predictor with the scan-independent part of the
    # current geometry
    self._reflection_predictor.update()

    # duck-typing for VaryingCrystalPredictionParameterisation. Only this
    # class has a compose(reflections) method. Sets ub_matrix (and caches
    # derivatives).
    try:
      self._prediction_parameterisation.compose(
          reflections)
    except AttributeError:
      pass

    # do prediction (updates reflection table in situ). Scan-varying prediction
    # is done automatically if the crystal has scan-points (assuming reflections
    # have ub_matrix set)

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

    # update xyzcal.px with the correct z_px values in keeping with above
    experiments = self._reflection_predictor._experiments
    for i, expt in enumerate(experiments):
      scan = expt.scan
      sel = (reflections['id'] == i)
      x_px, y_px, z_px = reflections['xyzcal.px'].select(sel).parts()
      if scan is not None:
        z_px = scan.get_array_index_from_angle(phi_calc.select(sel), deg=False)
      else:
        # must be a still image, z centroid not meaningful
        z_px = phi_calc.select(sel)
      xyzcal_px = flex.vec3_double(x_px, y_px, z_px)
      reflections['xyzcal.px'].set_selected(sel, xyzcal_px)

    # calculate residuals and assign columns
    reflections['x_resid'] = x_calc - x_obs
    reflections['x_resid2'] = reflections['x_resid']**2
    reflections['y_resid'] = y_calc - y_obs
    reflections['y_resid2'] = reflections['y_resid']**2
    reflections['phi_resid'] = phi_calc - phi_obs
    reflections['phi_resid2'] = reflections['phi_resid']**2

    return reflections

  def predict(self):
    """perform reflection prediction for the working reflections and update the
    reflection manager"""

    # get the matches
    reflections = self._reflection_manager.get_obs()

    # reset the 'use' flag for all observations
    self._reflection_manager.reset_accepted_reflections()

    # predict
    reflections = self._predict_core(reflections)

    # set used_in_refinement flag to all those that had predictions
    mask = reflections.get_flags(reflections.flags.predicted)
    reflections.set_flags(mask, reflections.flags.used_in_refinement)

    # collect the matches
    self.update_matches(force=True)

    return

  def predict_for_free_reflections(self):
    """perform prediction for the reflections not used for refinement"""

    refs = self._reflection_manager.get_free_reflections()
    if len(refs) == 0: return refs

    return self._predict_core(refs)

  def predict_for_reflection_table(self, reflections):
    """perform prediction for all reflections in the supplied table"""

    # set the entering flags if this has not been done
    from dials.algorithms.refinement.reflection_manager import calculate_entering_flags
    if not reflections.has_key("entering"):
      reflections['entering'] = calculate_entering_flags(reflections, self._experiments)

    # can only predict for experiments that exist and within the scan range
    # any other reflections will be left unchanged
    inc = flex.size_t_range(len(reflections))
    to_keep = flex.bool(len(inc), False)

    for iexp, exp in enumerate(self._experiments):
      sel = reflections['id'] == iexp

      # keep all reflections if there is no rotation axis
      if exp.goniometer is None:
        to_keep.set_selected(sel, True)
        continue

      # trim reflections outside the scan range
      phi = reflections['xyzobs.mm.value'].parts()[2]
      phi_min, phi_max = exp.scan.get_oscillation_range(deg=False)
      passed = (phi >= phi_min) & (phi <= phi_max)
      to_keep.set_selected(sel, passed)

    # determine indices to include and predict on the subset
    inc = inc.select(to_keep)
    sub_refl = reflections.select(inc)
    preds = self._predict_core(sub_refl)

    # set updated subset back into place
    reflections.set_selected(inc, preds)

    return reflections

  def calculate_gradients(self, reflections=None, callback=None):
    """delegate to the prediction_parameterisation object to calculate
    gradients for all the matched reflections, or just for those specified"""

    self.update_matches()

    if reflections:
      gradients = self._prediction_parameterisation.get_gradients(
        reflections, callback)
    else:
      gradients = self._prediction_parameterisation.get_gradients(
        self._matches, callback)

    return gradients

  def get_num_matches(self):
    """return the number of reflections currently used in the calculation"""

    self.update_matches()
    return len(self._matches)

  def get_num_matches_for_experiment(self, iexp=0):
    """return the number of reflections currently used in the calculation"""

    self.update_matches()
    sel = self._matches['id'] == iexp
    return sel.count(True)

  def get_num_matches_for_panel(self, ipanel=0):
    """return the number of reflections currently used in the calculation"""

    self.update_matches()
    sel = self._matches['panel'] == ipanel
    return sel.count(True)

  def update_matches(self, force=False):
    """ensure the observations matched to predictions are up to date"""

    if not self._matches or force:
      self._matches = self._reflection_manager.get_matches()

    return

  def compute_functional_gradients_and_curvatures(self, block=None):
    """calculate the value of the target function and its gradients. Set
    approximate curvatures as a side-effect"""

    self.update_matches()
    if block is not None:
      matches = block
    else:
      matches = self._matches
    nref = len(matches)

    # This is a hack for the case where nref=0. This should not be necessary
    # if bounds are provided for parameters to stop the algorithm exploring
    # unreasonable regions of parameter space where no predictions exist.
    # Unfortunately the L-BFGS line search does make such extreme trials.
    if nref == 0:
      self._curv = [1.] * len(self._prediction_parameterisation)
      return 1.e12, [1.] * len(self._prediction_parameterisation)

    residuals, weights = self._extract_residuals_and_weights(matches)
    w_resid = weights * residuals
    residuals2 = self._extract_squared_residuals(matches)

    # calculate target function
    L = 0.5 * flex.sum(weights * residuals2)

    def process_one_gradient(result):
      # copy gradients out of the result
      dX = result[self._grad_names[0]]
      dY = result[self._grad_names[1]]
      dZ = result[self._grad_names[2]]
      # reset result
      for k in result.keys():
        result[k] = None
      # add new keys
      grads = self._concatenate_gradients(dX, dY, dZ)
      result['dL_dp'] = flex.sum(w_resid * grads)
      result['curvature'] = flex.sum(weights * grads * grads)
      return result

    results = self.calculate_gradients(matches,
        callback=process_one_gradient)

    dL_dp = [result['dL_dp'] for result in results]
    curvs = [result['curvature'] for result in results]

    return (L, dL_dp, curvs)

  #def curvatures(self):
  #  """First order approximation to the diagonal of the Hessian based on the
  #  least squares form of the target"""
  #
  #  # relies on compute_functional_and_gradients being called first to set
  #  # self._curv
  #
  #  # Curvatures of zero will cause a crash, because their inverse is taken.
  #  assert all([c > 0.0 for c in self._curv])
  #
  #  return self._curv

  def compute_residuals(self):
    """return the vector of residuals plus their weights"""

    self.update_matches()
    return self._extract_residuals_and_weights(self._matches)

  def split_matches_into_blocks(self, nproc=1):
    """Return a list of the matches, split into blocks according to the
    gradient_calculation_blocksize parameter and the number of processes (if relevant).
    The number of blocks will be set such that the total number of reflections
    being processed by concurrent processes does not exceed gradient_calculation_blocksize"""

    self.update_matches()
    if self._gradient_calculation_blocksize:
      nblocks = int(ceil(len(self._matches) * nproc / self._gradient_calculation_blocksize))
    else:
      nblocks = nproc
    blocksize = int(ceil(len(self._matches) / nblocks))
    blocks = []
    for block_num in range(nblocks - 1):
      start = block_num * blocksize
      end = (block_num + 1) * blocksize
      blocks.append(self._matches[start:end])
    start = (nblocks - 1) * blocksize
    end = len(self._matches)
    blocks.append(self._matches[start:end])
    return blocks

  def compute_residuals_and_gradients(self, block=None):
    """return the vector of residuals plus their gradients and weights for
    non-linear least squares methods"""

    self.update_matches()
    if block is not None:
      matches = block
    else:
      matches = self._matches

    # Here we hardcode *three* types of residual, which might correspond to
    # X, Y, Phi (scans case) or X, Y, DeltaPsi (stills case).
    gradients = self.calculate_gradients(matches)
    dX_dp = [g[self._grad_names[0]] for g in gradients]
    dY_dp = [g[self._grad_names[1]] for g in gradients]
    dZ_dp = [g[self._grad_names[2]] for g in gradients]

    residuals, weights = self._extract_residuals_and_weights(matches)

    nelem = len(matches) * 3
    nparam = len(self._prediction_parameterisation)
    jacobian = self._build_jacobian(dX_dp, dY_dp, dZ_dp, nelem, nparam)

    return(residuals, jacobian, weights)

  def compute_restraints_residuals_and_gradients(self):
    '''delegate to the restraints_parameterisation object, if present, to
    calculate the vector of restraints residuals plus their gradients and
    weights for non-linear least squares methods'''

    if self._restraints_parameterisation:
      residuals, jacobian = \
        self._restraints_parameterisation.get_values_and_gradients()
      # currently weights are incorporated into the restraints, so provide
      # unit weights here
      weights = flex.double(len(residuals), 1.0)
      return(residuals, jacobian, weights)

    else:
      return None

  @staticmethod
  def _build_jacobian(dX_dp, dY_dp, dZ_dp, nelem, nparam):
    """construct Jacobian from lists of gradient vectors. This method may be
    overridden for the case where these vectors use sparse storage"""

    jacobian = flex.double(flex.grid(nelem, nparam))
    # loop over parameters
    for i in range(nparam):
      dX, dY, dZ = dX_dp[i], dY_dp[i], dZ_dp[i]
      col = flex.double.concatenate(dX, dY)
      col.extend(dZ)
      jacobian.matrix_paste_column_in_place(col, i)

    return jacobian

  @staticmethod
  def _concatenate_gradients(dX, dY, dZ):
    """concatenate three gradient vectors and return a flex.double. This method
    may be overriden for the case where these vectors use sparse storage"""

    grads = flex.double.concatenate(dX, dY)
    grads.extend(dZ)
    return grads

  @abc.abstractmethod
  def _extract_residuals_and_weights(matches):
    """extract vector of residuals and corresponding weights. The space the
    residuals are measured in (e.g. X, Y and Phi) and the order they are
    returned is determined by a concrete implementation of a staticmethod"""

    pass

  @abc.abstractmethod
  def _extract_squared_residuals(matches):
    """extract vector of squared residuals. The space the residuals are measured
    in (e.g. X, Y and Phi) and the order they are returned is determined by a
    concrete implementation of a staticmethod"""

    pass

  def rmsds(self):
    """calculate unweighted RMSDs for the matches"""

    self.update_matches()

    # cache rmsd calculation for achieved test
    self._rmsds = self._rmsds_core(self._matches)

    return self._rmsds

  def rmsds_for_reflection_table(self, reflections):
    """calculate unweighted RMSDs for the specified reflections. Caution: this
    assumes that the table reflections has the keys expected by _rmsds_core"""

    n = len(reflections)
    if n == 0: return None
    return self._rmsds_core(reflections)

  def rmsds_for_experiment(self, iexp=0):
    """calculate unweighted RMSDs for the selected experiment."""

    self.update_matches()
    sel = self._matches['id'] == iexp
    n = sel.count(True)
    if n == 0: return None

    rmsds = self._rmsds_core(self._matches.select(sel))
    return rmsds

  def rmsds_for_panel(self, ipanel=0):
    """calculate unweighted RMSDs for the selected panel."""

    self.update_matches()
    sel = self._matches['panel'] == ipanel
    n = sel.count(True)
    if n == 0: return None

    rmsds = self._rmsds_core(self._matches.select(sel))
    return rmsds

  @abc.abstractmethod
  def achieved(self):
    """return True to terminate the refinement."""

    pass

class LeastSquaresPositionalResidualWithRmsdCutoff(Target):
  """An implementation of the target class providing a least squares residual
  in terms of detector impact position X, Y and phi, terminating on achieved
  rmsd (or on intrisic convergence of the chosen minimiser)"""

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
    if experiments[0].scan:
      image_width_rad = abs(experiments[0].scan.get_oscillation(deg=False)[1])
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

  @staticmethod
  def _extract_residuals_and_weights(matches):

    # return residuals and weights as 1d flex.double vectors
    residuals = flex.double.concatenate(matches['x_resid'],
                                        matches['y_resid'])
    residuals.extend(matches['phi_resid'])

    weights, w_y, w_z = matches['xyzobs.mm.weights'].parts()
    weights.extend(w_y)
    weights.extend(w_z)

    return residuals, weights

  @staticmethod
  def _extract_squared_residuals(matches):

    residuals2 = flex.double.concatenate(matches['x_resid2'],
                                         matches['y_resid2'])
    residuals2.extend(matches['phi_resid2'])

    return residuals2

  def _rmsds_core(self, reflections):
    """calculate unweighted RMSDs for the specified reflections"""

    resid_x = flex.sum(reflections['x_resid2'])
    resid_y = flex.sum(reflections['y_resid2'])
    resid_phi = flex.sum(reflections['phi_resid2'])
    n = len(reflections)

    rmsds = (sqrt(resid_x / n),
             sqrt(resid_y / n),
             sqrt(resid_phi / n))
    return rmsds

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

class SparseGradientsMixin:
  """Mixin class to build a sparse Jacobian from gradients of the prediction
  formula stored as sparse vectors, and allow concatenation of gradient vectors
  that employed sparse storage."""

  @staticmethod
  def _build_jacobian(dX_dp, dY_dp, dZ_dp, nelem, nparam):
    """construct Jacobian from lists of sparse gradient vectors."""

    from scitbx import sparse
    nref = int(nelem / 3)
    X_mat = sparse.matrix(nref, nparam)
    Y_mat = sparse.matrix(nref, nparam)
    Z_mat = sparse.matrix(nref, nparam)
    jacobian = sparse.matrix(nelem, nparam)

    # loop over parameters, building full width blocks of the full Jacobian
    for i in range(nparam):
      X_mat[:,i] = dX_dp[i]
      Y_mat[:,i] = dY_dp[i]
      Z_mat[:,i] = dZ_dp[i]

    # set the blocks in the Jacobian
    jacobian.assign_block(X_mat, 0, 0)
    jacobian.assign_block(Y_mat, nref, 0)
    jacobian.assign_block(Z_mat, 2*nref, 0)

    return jacobian

  @staticmethod
  def _concatenate_gradients(dX, dY, dZ):
    """concatenate three sparse gradient vectors and return a flex.double."""

    dX = dX.as_dense_vector()
    dY = dY.as_dense_vector()
    dZ = dZ.as_dense_vector()
    grads = flex.double.concatenate(dX, dY)
    grads.extend(dZ)
    return grads

class LeastSquaresPositionalResidualWithRmsdCutoffSparse(
  SparseGradientsMixin, LeastSquaresPositionalResidualWithRmsdCutoff):
  """A version of the LeastSquaresPositionalResidualWithRmsdCutoff Target that
  uses a sparse matrix data structure for memory efficiency when there are a
  large number of Experiments"""

  pass
