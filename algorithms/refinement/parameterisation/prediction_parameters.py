#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
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

from dials.array_family import flex

"""This version of PredictionParameterisation vectorises calculations over
reflections, using flex arrays"""

class PredictionParameterisation(object):
  """
  Abstract interface for a class that groups together model parameterisations
  relating to diffraction geometry and provides:

  * A list of all free parameters concatenated from each of the models, with a
    getter and setter method that delegates to the contained models
  * Derivatives of the reflection prediction equation with respect to each of
    these free parameters

  Derived classes determine whether the reflection prediction equation is
  expressed in detector space (X, Y, phi) or orthogonalised reciprocal space.

  It is assumed that the provided model parameterisations will be one of four
  types:

  * Detector parameterisation
  * Beam parameterisation
  * Crystal orientation parameterisation
  * Crystal unit cell parameterisation

  One of each must be supplied for each Experiment. These could be satisfied by
  a dummy class if no parameterisation is desired for some model.

  We also need access to the underlying models that are parameterised. The
  model parameterisation objects do not provide access to these models as it is
  not their job to do so. Instead we keep a separate reference to an
  ExperimentList that allows access to the relevant models.

  The goniometer is not yet parameterised, but we need it for the equations if
  we are doing parameterisation in X, Y, Phi space. Conversely, if
  parameterisation is only in X, Y space, the goniometer model is optional.

  A class implementing PredictionParameterisation is used by a Refinery
  object directly, which takes the list of parameters, and indirectly via a
  Target function object, which takes the list of derivatives and composes the
  derivatives of a Target function from them.

  """

  def __init__(self,
               experiments,
               detector_parameterisations = None,
               beam_parameterisations = None,
               xl_orientation_parameterisations = None,
               xl_unit_cell_parameterisations = None):

    if detector_parameterisations is None:
      detector_parameterisations = []
    if beam_parameterisations is None:
      beam_parameterisations = []
    if xl_orientation_parameterisations is None:
      xl_orientation_parameterisations = []
    if xl_unit_cell_parameterisations is None:
      xl_unit_cell_parameterisations = []

    # References to the underlying models
    self._experiments = experiments

    # Keep references to all parameterised models
    self._detector_parameterisations = detector_parameterisations
    self._beam_parameterisations = beam_parameterisations
    self._xl_orientation_parameterisations = \
        xl_orientation_parameterisations
    self._xl_unit_cell_parameterisations = \
        xl_unit_cell_parameterisations

    # Check there are free parameters to refine
    self._length = self._len()
    if self._length == 0:
      raise RuntimeError("There are no free parameters for refinement")

    # Calculate Experiment to parameterisation mapping
    e2bp = {ids: i for i, dp in enumerate(beam_parameterisations) \
                 for ids in dp.get_experiment_ids()}
    e2xop = {ids: i for i, dp in enumerate(xl_orientation_parameterisations) \
                 for ids in dp.get_experiment_ids()}
    e2xucp = {ids: i for i, dp in enumerate(xl_unit_cell_parameterisations) \
                  for ids in dp.get_experiment_ids()}
    e2dp = {ids: i for i, dp in enumerate(detector_parameterisations) \
                 for ids in dp.get_experiment_ids()}
    from collections import namedtuple
    ParamSet = namedtuple('ParamSet', ['beam_param', 'xl_ori_param',
                                         'xl_uc_param', 'det_param'])

    self._exp_to_param = {i: ParamSet(e2bp.get(i), e2xop.get(i),
        e2xucp.get(i), e2dp.get(i)) for i, _ in enumerate(experiments)}

  def _len(self):

    length = 0
    for model in self._detector_parameterisations:
      length += model.num_free()
    for model in self._beam_parameterisations:
      length += model.num_free()
    for model in self._xl_orientation_parameterisations:
      length += model.num_free()
    for model in self._xl_unit_cell_parameterisations:
      length += model.num_free()

    return length

  def __len__(self):
    return self._length

  def get_param_vals(self):
    """return a concatenated list of parameters from each of the components
    in the global model"""

    global_p_list = []
    if self._detector_parameterisations:
      det_plists = [x.get_param_vals() for x
                    in self._detector_parameterisations]
      params = [x for l in det_plists for x in l]
      global_p_list.extend(params)

    if self._beam_parameterisations:
      src_plists = [x.get_param_vals() for x
                    in self._beam_parameterisations]
      params = [x for l in src_plists for x in l]
      global_p_list.extend(params)

    if self._xl_orientation_parameterisations:
      xlo_plists = [x.get_param_vals() for x
                    in self._xl_orientation_parameterisations]
      params = [x for l in xlo_plists for x in l]
      global_p_list.extend(params)

    if self._xl_unit_cell_parameterisations:
      xluc_plists = [x.get_param_vals() for x
                     in self._xl_unit_cell_parameterisations]
      params = [x for l in xluc_plists for x in l]
      global_p_list.extend(params)

    return global_p_list

  def get_param_names(self):
    """Return a list of the names of parameters in the order they are
    concatenated. Useful for output to log files and debugging. Use 1-based
    indexing for indices in the names"""
    param_names = []
    if self._detector_parameterisations:
      det_param_name_lists = [x.get_param_names() for x in \
                         self._detector_parameterisations]
      names = ["Detector%d" % (i + 1) + x for i, l \
               in enumerate(det_param_name_lists) for x in l]
      param_names.extend(names)

    if self._beam_parameterisations:
      beam_param_name_lists = [x.get_param_names() for x in \
                         self._beam_parameterisations]
      params = ["Beam%d" % (i + 1) + x for i, l \
                in enumerate(beam_param_name_lists) for x in l]
      param_names.extend(params)

    if self._xl_orientation_parameterisations:
      xlo_param_name_lists = [x.get_param_names() for x
                    in self._xl_orientation_parameterisations]
      params = ["Crystal%d" % (i + 1) + x for i, l \
                in enumerate(xlo_param_name_lists) for x in l]
      param_names.extend(params)

    if self._xl_unit_cell_parameterisations:
      xluc_param_name_lists = [x.get_param_names() for x
                     in self._xl_unit_cell_parameterisations]
      params = ["Crystal%d" % (i + 1) + x for i, l \
                in enumerate(xluc_param_name_lists) for x in l]
      param_names.extend(params)

    return param_names

  def set_param_vals(self, vals):
    """Set the parameter values of the contained models to the values in
    vals. This list must be of the same length as the result of
    get_param_vals and must contain the parameter values in the same order!
    This order is to be maintained by any sensible refinement engine."""

    assert len(vals) == len(self)
    it = iter(vals)

    for model in self._detector_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_vals(tmp)

    for model in self._beam_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_vals(tmp)

    for model in self._xl_orientation_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_vals(tmp)

    for model in self._xl_unit_cell_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_vals(tmp)

  def set_param_esds(self, esds):
    """Set the estimated standard deviations of parameter values of the
    contained models to the values in esds. This list must be of the same length
    as the result of get_param_vals and must contain the parameter values in the
    same order! This order is to be maintained by any sensible refinement
    engine."""

    assert len(esds) == len(self)
    it = iter(esds)

    for model in self._detector_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_esds(tmp)

    for model in self._beam_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_esds(tmp)

    for model in self._xl_orientation_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_esds(tmp)

    for model in self._xl_unit_cell_parameterisations:
      tmp = [it.next() for i in range(model.num_free())]
      model.set_param_esds(tmp)

  def calculate_model_state_uncertainties(self, var_cov):
    """
    Take the variance-covariance matrix of all free parameters calculated by
    the minimisation engine. For each parameterisation in the global model,
    extract the subset of this matrix for the associated block of parameters.
    Pass this on to the relevant model parameterisation to calculate its own
    uncertainty of state."""

    i = 0
    for model in self._detector_parameterisations:
      n = model.num_free()
      sub = var_cov.matrix_copy_block(i, i, n, n)
      state_covs = model.calculate_state_uncertainties(sub)
      if state_covs is None: continue
      if len(state_covs) == 1:
        model.set_state_uncertainties(state_covs[0])
      else:
        for i_state, state_cov in enumerate(state_covs):
          model.set_state_uncertainties(state_cov, multi_state_elt=i_state)
      i += n

    for model in self._beam_parameterisations:
      n = model.num_free()
      sub = var_cov.matrix_copy_block(i, i, n, n)
      state_covs = model.calculate_state_uncertainties(sub)
      if state_covs is None: continue
      if len(state_covs) == 1:
        model.set_state_uncertainties(state_covs[0])
      else:
        for i_state, state_cov in enumerate(state_covs):
          model.set_state_uncertainties(state_cov, multi_state_elt=i_state)
      i += n

    for model in self._xl_orientation_parameterisations:
      n = model.num_free()
      sub = var_cov.matrix_copy_block(i, i, n, n)
      state_covs = model.calculate_state_uncertainties(sub)
      if state_covs is None: continue
      if len(state_covs) == 1:
        model.set_state_uncertainties(state_covs[0])
      else:
        for i_state, state_cov in enumerate(state_covs):
          model.set_state_uncertainties(state_cov, multi_state_elt=i_state)
      i += n

    for model in self._xl_unit_cell_parameterisations:
      n = model.num_free()
      sub = var_cov.matrix_copy_block(i, i, n, n)
      state_covs = model.calculate_state_uncertainties(sub)
      if state_covs is None: continue
      if len(state_covs) == 1:
        model.set_state_uncertainties(state_covs[0])
      else:
        for i_state, state_cov in enumerate(state_covs):
          model.set_state_uncertainties(state_cov, multi_state_elt=i_state)
      i += n

    return

  def get_gradients(self, reflections, callback=None):
    """
    Calculate gradients of the prediction formula with respect to each
    of the parameters of the contained models, for all of the reflections.

    To be implemented by a derived class, which determines the space of the
    prediction formula (e.g. we calculate dX/dp, dY/dp, dphi/dp for the
    prediction formula for a rotation scan expressed in detector space, but
    components of d\vec{r}/dp for the prediction formula in reciprocal space

    """

    ### Calculate various quantities of interest for the reflections

    # Set up arrays of values for each reflection
    n = len(reflections)
    D = flex.mat3_double(n)
    s0 = flex.vec3_double(n)
    U = flex.mat3_double(n)
    B = flex.mat3_double(n)
    axis = flex.vec3_double(n)
    fixed_rotation = flex.mat3_double(n)

    for iexp, exp in enumerate(self._experiments):

      sel = reflections['id'] == iexp
      subref = reflections.select(sel)
      states = self._get_model_data_for_experiment(exp, subref)

      D.set_selected(sel, states['D'])
      s0.set_selected(sel, states['s0'])
      U.set_selected(sel, states['U'])
      B.set_selected(sel, states['B'])
      if exp.goniometer:
        axis.set_selected(sel, exp.goniometer.get_rotation_axis())
        fixed_rotation.set_selected(sel, exp.goniometer.get_fixed_rotation())

    return self._get_gradients_core(reflections, D, s0, U, B, axis, fixed_rotation, callback)

  @staticmethod
  def _extend_gradient_vectors(results, m, n,
                               keys=("dX_dp", "dY_dp", "dZ_dp")):
    """Extend results list by n empty results. These will each be a dictionary
    indexed by the given keys. The value for each key will be an empty vector of
    size m, to store the derivatives of n parameters, for m reflections.
    This method may be overriden by a derived class to e.g. use sparse
    vectors"""

    new_results = []
    for i in range(n):
      result = {}
      for key in keys:
        result[key] = flex.double(m, 0.)
      new_results.append(result)
    results.extend(new_results)

    return results

  def _get_model_data_for_experiment(self, experiment, reflections):
    """helper function to return model data s0, U, B and D for a particular
    experiment. D is always returned as an array the same length as the
    reflections for the experiment, whereas here U, B and s0 are returned as
    single matrices or vectors. In the scan-varying overload these will all be
    arrays."""

    # D matrix array
    D = flex.mat3_double(len(reflections))
    panels = reflections['panel']
    for ipanel, D_mat in enumerate([p.get_D_matrix() for p in experiment.detector]):
      sel = panels == ipanel
      D.set_selected(sel, D_mat)

    return {'s0':experiment.beam.get_s0(),
            'U':experiment.crystal.get_U(),
            'B':experiment.crystal.get_B(),
            'D':D}

  # The detector derivatives calculation is shared by scans and stills type
  # prediction, so this method is here, in the base class.
  def _detector_derivatives(self, pv, D, panel_id, parameterisation=None, dd_ddet_p=None):
    """helper function to convert derivatives of the detector state to
    derivatives of the vector pv. Derivatives that would all be null vectors
    are replaced with None"""

    if dd_ddet_p is None:

      # get the derivatives of detector d matrix for this panel
      dd_ddet_p = parameterisation.get_ds_dp(multi_state_elt=panel_id,
                                             use_none_as_null=True)

      # replace explicit null derivatives with None
      dd_ddet_p = [None if e is None else \
                   flex.mat3_double(len(D), e.elems) for e in dd_ddet_p]

    # calculate the derivative of pv for this parameter
    dpv_ddet_p = [der if der is None else (D * (der * -1.)) * pv for der in dd_ddet_p]

    return dpv_ddet_p

class SparseGradientVectorMixin(object):
  """Mixin class to use sparse vectors for storage of gradients of the
  prediction formula"""

  @staticmethod
  def _extend_gradient_vectors(results, m, n,
                               keys=("dX_dp", "dY_dp", "dZ_dp")):
    """Extend results list by n empty results. These will each be a dictionary
    indexed by the given keys. The value for each key will be an empty vector of
    size m, to store the derivatives of n parameters, for m reflections.
    This is the sparse vector version."""

    from scitbx import sparse

    new_results = []
    for i in range(n):
      result = {}
      for key in keys:
        result[key] = sparse.matrix_column(m)
      new_results.append(result)
    results.extend(new_results)

    return results

class XYPhiPredictionParameterisation(PredictionParameterisation):

  _grad_names = ("dX_dp", "dY_dp", "dphi_dp")

  def _get_gradients_core(self, reflections, D, s0, U, B, axis, fixed_rotation, callback=None):
    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    that reflects at rotation angle phi with scattering vector s that
    intersects panel panel_id. That is, calculate dX/dp, dY/dp and
    dphi/dp"""

    # Spindle rotation matrices for every reflection
    #R = self._axis.axis_and_angle_as_r3_rotation_matrix(phi)
    #R = flex.mat3_double(len(reflections))
    # NB for now use flex.vec3_double.rotate_around_origin each time I need the
    # rotation matrix R.

    self._axis = axis
    self._fixed_rotation = fixed_rotation
    self._s0 = s0

    # pv is the 'projection vector' for the ray along s1.
    self._D = D
    self._s1 = reflections['s1']
    self._pv = D * self._s1

    # also need quantities derived from pv, precalculated for efficiency
    u, v, w = self._pv.parts()
    self._w_inv = 1/w
    self._u_w_inv = u * self._w_inv
    self._v_w_inv = v * self._w_inv

    self._UB = U * B
    self._U = U
    self._B = B

    # r is the reciprocal lattice vector, in the lab frame
    self._h = reflections['miller_index'].as_vec3_double()
    self._phi_calc = reflections['xyzcal.mm'].parts()[2]
    self._r = (self._fixed_rotation * (self._UB * self._h)).rotate_around_origin(self._axis, self._phi_calc)

    # All of the derivatives of phi have a common denominator, given by
    # (e X r).s0, where e is the rotation axis. Calculate this once, here.
    self._e_X_r = self._axis.cross(self._r)
    self._e_r_s0 = (self._e_X_r).dot(self._s0)

    # Note that e_r_s0 -> 0 when the rotation axis, beam vector and
    # relp are coplanar. This occurs when a reflection just touches
    # the Ewald sphere.
    #
    # There is a relationship between e_r_s0 and zeta_factor.
    # Uncommenting the code below shows that
    # s0.(e X r) = zeta * |s X s0|

    #from dials.algorithms.profile_model.gaussian_rs import zeta_factor
    #from libtbx.test_utils import approx_equal
    #s = matrix.col(reflections['s1'][0])
    #z = zeta_factor(axis[0], s0[0], s)
    #ss0 = (s.cross(matrix.col(s0[0]))).length()
    #assert approx_equal(e_r_s0[0], z * ss0)

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

    # Set up empty list in which to store gradients
    m = len(reflections)
    results = []

    # determine experiment to indices mappings once, here
    experiment_to_idx = []
    for iexp, exp in enumerate(self._experiments):

      sel = reflections['id'] == iexp
      isel = sel.iselection()
      experiment_to_idx.append(isel)

    # reset a pointer to the parameter number
    self._iparam = 0

    ### Work through the parameterisations, calculating their contributions
    ### to derivatives d[pv]/dp and d[phi]/dp

    # loop over the detector parameterisations
    for dp in self._detector_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in dp.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Access the detector model being parameterised
      detector = dp.get_model()

      # Get panel numbers of the affected reflections
      panel = reflections['panel'].select(isel)

      # Extend derivative vectors for this detector parameterisation
      results = self._extend_gradient_vectors(results, m, dp.num_free(),
        keys=self._grad_names)

      # loop through the panels in this detector
      for panel_id, _ in enumerate(exp.detector):

        # get the right subset of array indices to set for this panel
        sub_isel = isel.select(panel == panel_id)
        if len(sub_isel) == 0:
          # if no reflections intersect this panel, skip calculation
          continue
        sub_pv = self._pv.select(sub_isel)
        sub_D = self._D.select(sub_isel)
        dpv_ddet_p = self._detector_derivatives(sub_pv, sub_D, panel_id,
          parameterisation=dp)

        # convert to dX/dp, dY/dp and assign the elements of the vectors
        # corresponding to this experiment and panel
        sub_w_inv = self._w_inv.select(sub_isel)
        sub_u_w_inv = self._u_w_inv.select(sub_isel)
        sub_v_w_inv = self._v_w_inv.select(sub_isel)
        dX_ddet_p, dY_ddet_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
          sub_w_inv, sub_u_w_inv, sub_v_w_inv, dpv_ddet_p)

        # use a local parameter index pointer because we set all derivatives
        # for this panel before moving on to the next
        iparam = self._iparam
        for dX, dY in zip(dX_ddet_p, dY_ddet_p):
          if dX is not None:
            results[iparam][self._grad_names[0]].set_selected(sub_isel, dX)
          if dY is not None:
            results[iparam][self._grad_names[1]].set_selected(sub_isel, dY)
          # increment the local parameter index pointer
          iparam += 1

      if callback is not None:
        iparam = self._iparam
        for i in range(dp.num_free()):
          results[iparam] = callback(results[iparam])
          iparam += 1

      # increment the parameter index pointer to the last detector parameter
      self._iparam += dp.num_free()

    # loop over the beam parameterisations
    for bp in self._beam_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in bp.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Extend derivative vectors for this beam parameterisation
      results = self._extend_gradient_vectors(results, m, bp.num_free(),
        keys=self._grad_names)

      if len(isel) == 0:
        # if no reflections are in this experiment, skip calculation of
        # gradients, but must still process null gradients by a callback
        if callback is not None:
          for iparam in xrange(bp.num_free()):
            results[self._iparam] = callback(results[self._iparam])
            self._iparam += 1
        else:
          self._iparam += bp.num_free()
        continue

      # Get required data from those reflections
      r = self._r.select(isel)
      e_X_r = self._e_X_r.select(isel)
      e_r_s0 = self._e_r_s0.select(isel)
      D = self._D.select(isel)

      w_inv = self._w_inv.select(isel)
      u_w_inv = self._u_w_inv.select(isel)
      v_w_inv = self._v_w_inv.select(isel)

      dpv_dbeam_p, dphi_dbeam_p = self._beam_derivatives(r, e_X_r, e_r_s0, D,
        parameterisation=bp)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dbeam_p, dY_dbeam_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dbeam_p)
      for dX, dY, dphi in zip(dX_dbeam_p, dY_dbeam_p, dphi_dbeam_p):
        results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        results[self._iparam][self._grad_names[2]].set_selected(isel, dphi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    # loop over the crystal orientation parameterisations
    for xlop in self._xl_orientation_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in xlop.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Extend derivative vectors for this crystal orientation parameterisation
      results = self._extend_gradient_vectors(results, m, xlop.num_free(),
        keys=self._grad_names)

      if len(isel) == 0:
        # if no reflections are in this experiment, skip calculation of
        # gradients, but must still process null gradients by a callback
        if callback is not None:
          for iparam in xrange(xlop.num_free()):
            results[self._iparam] = callback(results[self._iparam])
            self._iparam += 1
        else:
          self._iparam += xlop.num_free()
        continue

      # Get required data from those reflections
      axis = self._axis.select(isel)
      fixed_rotation = self._fixed_rotation.select(isel)
      phi_calc = self._phi_calc.select(isel)
      h = self._h.select(isel)
      s1 = self._s1.select(isel)
      e_X_r = self._e_X_r.select(isel)
      e_r_s0 = self._e_r_s0.select(isel)
      B = self._B.select(isel)
      D = self._D.select(isel)

      w_inv = self._w_inv.select(isel)
      u_w_inv = self._u_w_inv.select(isel)
      v_w_inv = self._v_w_inv.select(isel)

      dpv_dxlo_p, dphi_dxlo_p = self._xl_orientation_derivatives(
          axis, fixed_rotation, phi_calc, h, s1, e_X_r, e_r_s0, B, D,
          parameterisation=xlop)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dxlo_p, dY_dxlo_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dxlo_p)
      for dX, dY, dphi in zip(dX_dxlo_p, dY_dxlo_p, dphi_dxlo_p):
        if dX is not None:
          results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        if dY is not None:
          results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        if dphi is not None:
          results[self._iparam][self._grad_names[2]].set_selected(isel, dphi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    # loop over the crystal unit cell parameterisations
    for xlucp in self._xl_unit_cell_parameterisations:

      # Determine (sub)set of reflections affected by this parameterisation
      isel = flex.size_t()
      for exp_id in xlucp.get_experiment_ids():
        isel.extend(experiment_to_idx[exp_id])

      # Extend derivative vectors for this crystal unit cell parameterisation
      results = self._extend_gradient_vectors(results, m, xlucp.num_free(),
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

      # Get required data from those reflections
      axis = self._axis.select(isel)
      fixed_rotation = self._fixed_rotation.select(isel)
      phi_calc = self._phi_calc.select(isel)
      h = self._h.select(isel)
      s1 = self._s1.select(isel)
      e_X_r = self._e_X_r.select(isel)
      e_r_s0 = self._e_r_s0.select(isel)
      U = self._U.select(isel)
      D = self._D.select(isel)

      w_inv = self._w_inv.select(isel)
      u_w_inv = self._u_w_inv.select(isel)
      v_w_inv = self._v_w_inv.select(isel)

      dpv_dxluc_p, dphi_dxluc_p =  self._xl_unit_cell_derivatives(
        axis, fixed_rotation, phi_calc, h, s1, e_X_r, e_r_s0, U, D,
        parameterisation=xlucp)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dxluc_p, dY_dxluc_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dxluc_p)
      for dX, dY, dphi in zip(dX_dxluc_p, dY_dxluc_p, dphi_dxluc_p):
        if dX is not None:
          results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        if dY is not None:
          results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        if dphi is not None:
          results[self._iparam][self._grad_names[2]].set_selected(isel, dphi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    return results

  def _beam_derivatives(self, r, e_X_r, e_r_s0, D, parameterisation=None, ds0_dbeam_p=None):
    """helper function to extend the derivatives lists by
    derivatives of the beam parameterisations"""

    if ds0_dbeam_p is None:

      # get the derivatives of the beam vector wrt the parameters
      ds0_dbeam_p = parameterisation.get_ds_dp(use_none_as_null=True)

      ds0_dbeam_p = [None if e is None else flex.vec3_double(len(r), e.elems) \
                     for e in ds0_dbeam_p]

    dphi_dp = []
    dpv_dp = []

    # loop through the parameters
    for der in ds0_dbeam_p:

      if der is None:
        dphi_dp.append(None)
        dpv_dp.append(None)
        continue

      # calculate the derivative of phi for this parameter
      dphi = (r.dot(der) / e_r_s0) * -1.0
      dphi_dp.append(dphi)

      # calculate the derivative of pv for this parameter
      dpv_dp.append(D * (e_X_r * dphi + der))

    return dpv_dp, dphi_dp

  def _xl_orientation_derivatives(self, axis, fixed_rotation, phi_calc, h,
    s1, e_X_r, e_r_s0, B, D, parameterisation=None, dU_dxlo_p=None):
    """helper function to extend the derivatives lists by
    derivatives of the crystal orientation parameterisations"""

    if dU_dxlo_p is None:

      # get derivatives of the U matrix wrt the parameters
      dU_dxlo_p = [None if der is None else flex.mat3_double(len(B), der.elems) \
                   for der in parameterisation.get_ds_dp(use_none_as_null=True)]

    dphi_dp = []
    dpv_dp = []

    # loop through the parameters
    for der in dU_dxlo_p:

      if der is None:
        dphi_dp.append(None)
        dpv_dp.append(None)
        continue

      # calculate the derivative of r for this parameter
      # FIXME COULD DO THIS BETTER WITH __rmul__?!
      tmp = fixed_rotation * (der * B * h)
      dr = tmp.rotate_around_origin(axis, phi_calc)

      # calculate the derivative of phi for this parameter
      dphi = -1.0 * dr.dot(s1) / e_r_s0
      dphi_dp.append(dphi)

      # calculate the derivative of pv for this parameter
      dpv_dp.append(D * (dr + e_X_r * dphi))

    return dpv_dp, dphi_dp

  def _xl_unit_cell_derivatives(self, axis, fixed_rotation, phi_calc, h,
    s1, e_X_r, e_r_s0, U, D, parameterisation=None, dB_dxluc_p=None):
    """helper function to extend the derivatives lists by
    derivatives of the crystal unit cell parameterisations"""

    if dB_dxluc_p is None:

      # get derivatives of the B matrix wrt the parameters
      dB_dxluc_p = [None if der is None else flex.mat3_double(len(U), der.elems) \
                    for der in parameterisation.get_ds_dp(use_none_as_null=True)]

    dphi_dp = []
    dpv_dp = []

    # loop through the parameters
    for der in dB_dxluc_p:

      if der is None:
        dphi_dp.append(None)
        dpv_dp.append(None)
        continue

      # calculate the derivative of r for this parameter
      tmp = fixed_rotation * (U * der * h)
      dr = tmp.rotate_around_origin(axis, phi_calc)

      # calculate the derivative of phi for this parameter
      dphi = -1.0 * dr.dot(s1) / e_r_s0
      dphi_dp.append(dphi)

      # calculate the derivative of pv for this parameter
      dpv_dp.append(D * (dr + e_X_r * dphi))

    return dpv_dp, dphi_dp

  @staticmethod
  def _calc_dX_dp_and_dY_dp_from_dpv_dp(w_inv, u_w_inv, v_w_inv, dpv_dp):
    """helper function to calculate positional derivatives from
    dpv_dp using the quotient rule"""

    dX_dp = []
    dY_dp = []

    for der in dpv_dp:
      if der is None:
        dX_dp.append(None)
        dY_dp.append(None)
      else:
        du_dp, dv_dp, dw_dp = der.parts()

        dX_dp.append(w_inv * (du_dp - dw_dp * u_w_inv))
        dY_dp.append(w_inv * (dv_dp - dw_dp * v_w_inv))

    return dX_dp, dY_dp

class XYPhiPredictionParameterisationSparse(SparseGradientVectorMixin,
  XYPhiPredictionParameterisation):
  """A version of XYPhiPredictionParameterisation that uses a sparse matrix
  data structure for memory efficiency when there are a large number of
  Experiments"""
  pass
