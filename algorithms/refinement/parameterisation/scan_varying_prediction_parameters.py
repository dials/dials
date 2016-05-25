#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# Python imports
from __future__ import division
from math import floor

# cctbx imports
from scitbx import matrix

# dials imports
#from dials_refinement_helpers_ext import *
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation, SparseGradientVectorMixin

from dials.array_family import flex

class ScanVaryingPredictionParameterisation(XYPhiPredictionParameterisation):
  """Support model parameterisations that vary with time (via the proxy of
  "observed image number")"""

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

    # determine once which types of parameterisations are scan-varying
    self._varying_detectors = any(hasattr(p, 'num_sets') for p in detector_parameterisations)
    self._varying_beams = any(hasattr(p, 'num_sets') for p in beam_parameterisations)
    self._varying_xl_orientations = any(hasattr(p, 'num_sets') for p in xl_orientation_parameterisations)
    self._varying_xl_unit_cells = any(hasattr(p, 'num_sets') for p in xl_unit_cell_parameterisations)

    super(ScanVaryingPredictionParameterisation, self).__init__(
      experiments,
      detector_parameterisations = detector_parameterisations,
      beam_parameterisations = beam_parameterisations,
      xl_orientation_parameterisations = xl_orientation_parameterisations,
      xl_unit_cell_parameterisations = xl_unit_cell_parameterisations)

  def _get_xl_orientation_parameterisation(self, experiment_id):
    """Return the crystal orientation parameterisation for the requested
    experiment number (or None if the crystal orientation in that experiment
    is not parameterised)"""

    param_set = self._exp_to_param[experiment_id]
    xl_op = None
    if param_set.xl_ori_param is not None:
      xl_op = self._xl_orientation_parameterisations[param_set.xl_ori_param]

    return xl_op

  def _get_xl_unit_cell_parameterisation(self, experiment_id):
    """Return the crystal unit cell parameterisation for the requested
    experiment number (or None if the crystal unit cell in that experiment
    is not parameterised)"""

    param_set = self._exp_to_param[experiment_id]
    xl_ucp = None
    if param_set.xl_uc_param is not None:
      xl_ucp = self._xl_unit_cell_parameterisations[param_set.xl_uc_param]

    return xl_ucp

  def _get_beam_parameterisation(self, experiment_id):
    """Return the beam parameterisation for the requested
    experiment number (or None if the beam in that experiment
    is not parameterised)"""

    param_set = self._exp_to_param[experiment_id]
    bp = None
    if param_set.beam_param is not None:
      bp = self._beam_parameterisations[param_set.beam_param]

    return bp

  def _get_detector_parameterisation(self, experiment_id):
    """Return the detector parameterisation for the requested
    experiment number (or None if the detector in that experiment
    is not parameterised)"""

    param_set = self._exp_to_param[experiment_id]
    dp = None
    if param_set.det_param is not None:
      dp = self._detector_parameterisations[param_set.det_param]

    return dp

  def _get_state_from_parameterisation(self, parameterisation, frame,
                                       multi_state_elt=None):
    """Get the model state from the parameterisation at the specified frame,
    taking care of whether it is a scan-varying parameterisation or not"""

    if parameterisation is None: return None
    if hasattr(parameterisation, 'num_sets'):
      parameterisation.compose(frame)
    if multi_state_elt is None:
      state = parameterisation.get_state()
    else:
      state = parameterisation.get_state(multi_state_elt = multi_state_elt)
    return state

  def _prepare_for_compose(self, reflections):

    nref = len(reflections)
    # set columns if needed
    if 'u_matrix' not in reflections:
      reflections['u_matrix'] = flex.mat3_double(nref)
    if 'b_matrix' not in reflections:
      reflections['b_matrix'] = flex.mat3_double(nref)
    if 's0_vector' not in reflections:
      reflections['s0_vector'] = flex.vec3_double(nref)
    if 'd_matrix' not in reflections:
      reflections['d_matrix'] = flex.mat3_double(nref)

    # set columns in the reflection table to store the derivative of state for
    # each reflection, if needed
    null9 = (0., 0., 0., 0., 0., 0., 0., 0., 0.)
    null3 = (0., 0., 0.)
    if self._varying_xl_orientations and "dU_dp0" not in reflections:
      max_free_params = max([e.num_free() for e in self._xl_orientation_parameterisations])
      for i in range(max_free_params):
        colname = "dU_dp{0}".format(i)
        reflections[colname] = flex.mat3_double(nref, null9)
    if self._varying_xl_unit_cells and "dB_dp0" not in reflections:
      max_free_params = max([e.num_free() for e in self._xl_unit_cell_parameterisations])
      for i in range(max_free_params):
        colname = "dB_dp{0}".format(i)
        reflections[colname] = flex.mat3_double(nref, null9)
    if self._varying_detectors and "dd_dp0" not in reflections:
      max_free_params = max([e.num_free() for e in self._detector_parameterisations])
      for i in range(max_free_params):
        colname = "dd_dp{0}".format(i)
        reflections[colname] = flex.mat3_double(nref, null9)
    if self._varying_beams and "ds0_dp0" not in reflections:
      max_free_params = max([e.num_free() for e in self._beam_parameterisations])
      for i in range(max_free_params):
        colname = "ds0_dp{0}".format(i)
        reflections[colname] = flex.vec3_double(nref, null3)

    return

  def compose(self, reflections):
    """Compose scan-varying crystal parameterisations at the specified image
    number, for the specified experiment, for all reflections. Put the U, B and
    UB matrices in the reflection table, and cache the derivatives."""

    self._prepare_for_compose(reflections)

    for iexp, exp in enumerate(self._experiments):

      # select the reflections of interest
      sel = reflections['id'] == iexp
      isel = sel.iselection()

      # get their frame numbers and intersecting panels
      obs_image_numbers = (reflections['xyzobs.px.value'].parts()[2]).select(isel)
      panels = reflections['panel'].select(isel)

      # identify which parameterisations to use for this experiment
      xl_op = self._get_xl_orientation_parameterisation(iexp)
      xl_ucp = self._get_xl_unit_cell_parameterisation(iexp)
      bp = self._get_beam_parameterisation(iexp)
      dp = self._get_detector_parameterisation(iexp)

      # get state and derivatives for each reflection
      for i, frame, pnl in zip(isel, obs_image_numbers, panels):

        # model states at current frame
        U = self._get_state_from_parameterisation(xl_op, frame)
        if U is None: U = exp.crystal.get_U()

        B = self._get_state_from_parameterisation(xl_ucp, frame)
        if B is None: B = exp.crystal.get_B()

        s0 = self._get_state_from_parameterisation(bp, frame)
        if s0 is None: s0 = exp.beam.get_s0()

        if dp is not None and dp.is_multi_state():
          dmat  = self._get_state_from_parameterisation(dp, frame, multi_state_elt = panel)
        else:
          dmat  = self._get_state_from_parameterisation(dp, frame)
        if dmat is None: dmat = exp.detector[panel].get_d_matrix()

        # set states and their derivatives into reflections
        row = {'u_matrix':U.elems,
               'b_matrix':B.elems,
               's0_vector':s0,
               'd_matrix':dmat}
        if xl_op is not None and self._varying_xl_orientations:
          for j, dU in enumerate(xl_op.get_ds_dp()):
            colname = "dU_dp{0}".format(j)
            row[colname] = dU
        if xl_ucp is not None and self._varying_xl_unit_cells:
          for j, dB in enumerate(xl_ucp.get_ds_dp()):
            colname = "dB_dp{0}".format(j)
            row[colname] = dB
        if bp is not None and self._varying_beams:
          for j, ds0 in enumerate(bp.get_ds_dp()):
            colname = "ds0_dp{0}".format(j)
            row[colname] = ds0
        if dp is not None and self._varying_detectors:
          if dp.is_multi_state:
            dds = dp.get_ds_dp(multi_state_elt=pnl)
          else:
            dds = dp.get_ds_dp()
          for j, dd in enumerate(dds):
            colname = "dd_dp{0}".format(j)
            row[colname] = dd
        reflections[i] = row

    # set the UB matrices for prediction
    reflections['ub_matrix'] = reflections['u_matrix'] * reflections['b_matrix']

    return

  def get_UB(self, obs_image_number, experiment_id):
    """Extract the setting matrix from the contained scan-dependent crystal
    parameterisations at specified image number."""

    # called by refiner.run for setting the crystal scan points

    # identify which crystal parameterisations to use for this experiment
    xl_op = self._get_xl_orientation_parameterisation(experiment_id)
    xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

    # model states at current frame
    U = self._get_state_from_parameterisation(xl_op, obs_image_number)
    if U is None: U = self._experiments[experiment_id].crystal.get_U()

    B = self._get_state_from_parameterisation(xl_ucp, obs_image_number)
    if B is None: B = self._experiments[experiment_id].crystal.get_B()

    return U*B

  # overloaded for the scan-varying case
  def _get_U_B_for_experiment(self, crystal, reflections, isel):
    """helper function to return either a single U, B pair (for scan-static) or
    U, B arrays (scan-varying; overloaded in derived class) for a particular
    experiment."""

    # crystal ignored here (it is needed for the scan-static version only)
    U = reflections['u_matrix'].select(isel)
    B = reflections['b_matrix'].select(isel)
    return U, B

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
            results[iparam]['dX_dp'].set_selected(sub_isel, dX)
          if dY is not None:
            results[iparam]['dY_dp'].set_selected(sub_isel, dY)
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
        # if no reflections are in this experiment, skip calculation
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
        # if no reflections are in this experiment, skip calculation
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

      # get derivatives of the U matrix wrt the parameters
      dU_dxlo_p = [reflections["dU_dp{0}".format(i)].select(isel) \
                   for i in range(xlop.num_free())]
      dpv_dxlo_p, dphi_dxlo_p = self._xl_orientation_derivatives(
        axis, fixed_rotation, phi_calc, h, s1, e_X_r, e_r_s0, B, D,
        dU_dxlo_p=dU_dxlo_p)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dxlo_p, dY_dxlo_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dxlo_p)
      for dX, dY, dphi in zip(dX_dxlo_p, dY_dxlo_p, dphi_dxlo_p):
        results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
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
        # if no reflections are in this experiment, skip calculation
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

      dB_dxluc_p = [reflections["dB_dp{0}".format(i)].select(isel) \
                   for i in range(xlucp.num_free())]
      dpv_dxluc_p, dphi_dxluc_p =  self._xl_unit_cell_derivatives(
        axis, fixed_rotation, phi_calc, h, s1, e_X_r, e_r_s0, U, D,
        dB_dxluc_p=dB_dxluc_p)

      # convert to dX/dp, dY/dp and assign the elements of the vectors
      # corresponding to this experiment
      dX_dxluc_p, dY_dxluc_p = self._calc_dX_dp_and_dY_dp_from_dpv_dp(
        w_inv, u_w_inv, v_w_inv, dpv_dxluc_p)
      for dX, dY, dphi in zip(dX_dxluc_p, dY_dxluc_p, dphi_dxluc_p):
        results[self._iparam][self._grad_names[0]].set_selected(isel, dX)
        results[self._iparam][self._grad_names[1]].set_selected(isel, dY)
        results[self._iparam][self._grad_names[2]].set_selected(isel, dphi)
        if callback is not None:
          results[self._iparam] = callback(results[self._iparam])
        # increment the parameter index pointer
        self._iparam += 1

    return results

  def calculate_model_state_uncertainties(self, var_cov=None,
                                          obs_image_number=None,
                                          experiment_id=None):
    """
    Take the variance-covariance matrix of all free parameters calculated by
    the minimisation engine. For each parameterisation in the global model,
    extract the subset of this matrix for the associated block of parameters.
    Pass this on to the relevant model parameterisation to calculate its own
    uncertainty of state.

    This scan-varying version should first be called with var_cov set but
    obs_image_number=None and experiment_id=None. This calls the scan-static
    version to do the calculation for the scan-static parameterisations and
    also caches the subsets of var_cov relevant for the scan-varying
    parameterisations. Subsequent calls should provide obs_image_number and
    experiment_id to calculate for a particular crystal at a particular
    scan-point"""

    # first call, only a variance-covariance matrix is supplied
    if var_cov is not None:
      assert [obs_image_number, experiment_id].count(None) == 2
      super(ScanVaryingPredictionParameterisation,
            self).calculate_model_state_uncertainties(var_cov)
      return

    # later calls, only an experiment and image number are supplied
    else:
      # identify the crystal parameterisations for this experiment
      xl_op = self._get_xl_orientation_parameterisation(experiment_id)
      xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

      # compose at the requested image number and calculate using the cached
      # varcov matrices. Take the first elt of the list becase the crystal
      # parameterisations are not multi-state
      try:
        xl_op.compose(obs_image_number)
        U_cov = xl_op.calculate_state_uncertainties(var_cov=None)[0]
      except AttributeError:
        U_cov = None

      try:
        xl_ucp.compose(obs_image_number)
        B_cov = xl_ucp.calculate_state_uncertainties(var_cov=None)[0]
      except AttributeError:
        B_cov = None

    return U_cov, B_cov

  def set_model_state_uncertainties(self, u_cov_list, b_cov_list,
                                          experiment_id=None):
    """Identify the parameterisation"""

    xl_op = self._get_xl_orientation_parameterisation(experiment_id)
    xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)
    try:
      xl_op.set_state_uncertainties(u_cov_list)
    except AttributeError:
      pass

    try:
      xl_ucp.set_state_uncertainties(b_cov_list)
    except AttributeError:
      pass

    return

class ScanVaryingPredictionParameterisationFast(ScanVaryingPredictionParameterisation):
  """Overloads compose to calculate model states per block rather than per
  reflection"""

  def compose(self, reflections):
    """Compose scan-varying crystal parameterisations at the specified image
    number, for the specified experiment, for each image. Put the U, B and
    UB matrices in the reflection table, and cache the derivatives."""

    self._prepare_for_compose(reflections)

    for iexp, exp in enumerate(self._experiments):

      # select the reflections of interest
      sel = reflections['id'] == iexp
      isel = sel.iselection()

      blocks = reflections['block'].select(isel)

      # identify which parameterisations to use for this experiment
      xl_op = self._get_xl_orientation_parameterisation(iexp)
      xl_ucp = self._get_xl_unit_cell_parameterisation(iexp)
      bp = self._get_beam_parameterisation(iexp)
      dp = self._get_detector_parameterisation(iexp)

      # get state and derivatives for each block
      for block in xrange(flex.min(blocks),
                          flex.max(blocks) + 1):

        # determine the subset of reflections this affects
        subsel = isel.select(blocks == block)
        if len(subsel) == 0: continue

        # get the panels hit by these reflections
        panels = reflections['panel'].select(subsel)

        # get the integer frame number nearest the centre of that block
        frames = reflections['block_centre'].select(subsel)

        # can only be false if original block assignment has gone wrong
        assert frames.all_eq(frames[0]), \
            "Failing: a block contains reflections that shouldn't be there"
        frame = int(floor(frames[0]))

        # model states at current frame
        U = self._get_state_from_parameterisation(xl_op, frame)
        if U is None: U = exp.crystal.get_U()

        B = self._get_state_from_parameterisation(xl_ucp, frame)
        if B is None: B = exp.crystal.get_B()

        s0 = self._get_state_from_parameterisation(bp, frame)
        if s0 is None: s0 = exp.beam.get_s0()

        # set states for crystal and beam
        reflections['u_matrix'].set_selected(subsel, U.elems)
        reflections['b_matrix'].set_selected(subsel, B.elems)
        reflections['s0_vector'].set_selected(subsel, s0.elems)

        # set states and derivatives for multi-panel detector
        if dp is not None and dp.is_multi_state():

          # loop through the panels in this detector
          for panel_id, _ in enumerate(exp.detector):

            # get the right subset of array indices to set for this panel
            subsel2 = subsel.select(panels == panel_id)
            if len(subsel2) == 0:
              # if no reflections intersect this panel, skip calculation
              continue

            dmat  = self._get_state_from_parameterisation(dp, frame, multi_state_elt=panel_id)
            if dmat is None: dmat = exp.detector[panel_id].get_d_matrix()
            reflections['d_matrix'].set_selected(subsel2, dmat)

            if dp is not None and self._varying_detectors:
              for j, dd in enumerate(dp.get_ds_dp(multi_state_elt=panel_id)):
                colname = "dd_dp{0}".format(j)
                reflections[colname].set_selected(subsel, dd)

        else: # set states and derivatives for single panel detector

          dmat  = self._get_state_from_parameterisation(dp, frame)
          if dmat is None: dmat = exp.detector[panel].get_d_matrix()
          reflections['d_matrix'].set_selected(subsel, dmat)

          if dp is not None and self._varying_detectors:
            for j, dd in enumerate(dp.get_ds_dp()):
              colname = "dd_dp{0}".format(j)
              reflections[colname].set_selected(subsel, dd)

        # set derivatives of the states for crystal and beam
        if xl_op is not None and self._varying_xl_orientations:
          for j, dU in enumerate(xl_op.get_ds_dp()):
            colname = "dU_dp{0}".format(j)
            reflections[colname].set_selected(subsel, dU)
        if xl_ucp is not None and self._varying_xl_unit_cells:
          for j, dB in enumerate(xl_ucp.get_ds_dp()):
            colname = "dB_dp{0}".format(j)
            reflections[colname].set_selected(subsel, dB)
        if bp is not None and self._varying_beams:
          for j, ds0 in enumerate(bp.get_ds_dp()):
            colname = "ds0_dp{0}".format(j)
            reflections[colname].set_selected(subsel, ds0)

    # set the UB matrices for prediction
    reflections['ub_matrix'] = reflections['u_matrix'] * reflections['b_matrix']

    return

class ScanVaryingPredictionParameterisationSparse(
    SparseGradientVectorMixin, ScanVaryingPredictionParameterisation):
  pass

class ScanVaryingPredictionParameterisationFastSparse(
    SparseGradientVectorMixin, ScanVaryingPredictionParameterisationFast):
  pass
