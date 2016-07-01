#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from math import floor
from scitbx import matrix
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation, SparseGradientVectorMixin

from dials.array_family import flex

class ScanVaryingPredictionParameterisation(XYPhiPredictionParameterisation):
  """An extension of the rotation scans version of the
  PredictionParameterisation class that supports model parameterisations that
  vary smoothly with the observed image number"""

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

    # determine once here which types of parameterisations are scan-varying
    self._varying_detectors = any(hasattr(p, 'num_sets')
      for p in detector_parameterisations)
    self._varying_beams = any(hasattr(p, 'num_sets')
      for p in beam_parameterisations)
    self._varying_xl_orientations = any(hasattr(p, 'num_sets')
      for p in xl_orientation_parameterisations)
    self._varying_xl_unit_cells = any(hasattr(p, 'num_sets')
      for p in xl_unit_cell_parameterisations)

    # set up base class
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
    """Return the beam parameterisation for the requested experiment number
    (or None if the beam in that experiment is not parameterised)"""

    param_set = self._exp_to_param[experiment_id]
    bp = None
    if param_set.beam_param is not None:
      bp = self._beam_parameterisations[param_set.beam_param]

    return bp

  def _get_detector_parameterisation(self, experiment_id):
    """Return the detector parameterisation for the requested experiment number
    (or None if the detector in that experiment is not parameterised)"""

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
    if hasattr(parameterisation, 'num_sets') and not self._current_frame.get(
        parameterisation) == frame:
      parameterisation.compose(frame)
      self._current_frame[parameterisation] = frame
    if multi_state_elt is None:
      state = parameterisation.get_state()
    else:
      state = parameterisation.get_state(multi_state_elt = multi_state_elt)
    return state

  def _prepare_for_compose(self, reflections):
    """Add columns to the reflection table to hold the varying state matrices
    or vectors for the experimental models, if required. Also add columns for
    the derivatives of states that are scan-varying"""

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
    if 'D_matrix' not in reflections:
      reflections['D_matrix'] = flex.mat3_double(nref)

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

  def compose(self, reflections, skip_derivatives=False):
    """Compose scan-varying crystal parameterisations at each observed image
    number, for each experiment, for all reflections. Put the U, B and UB
    matrices in the reflection table, and cache the derivatives if requested"""

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

      # reset current frame cache for scan-varying parameterisations
      self._current_frame = {}

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
        Dmat = exp.detector[panel].get_D_matrix() # actually need D, but need
        # the above to compose at the right frame

        # set states and their derivatives into reflections
        row = {'u_matrix':U.elems,
               'b_matrix':B.elems,
               's0_vector':s0,
               'd_matrix':dmat,
               'D_matrix':Dmat}
        if not skip_derivatives:
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

  # called by refiner.run for setting the crystal scan points
  def get_UB(self, obs_image_number, experiment_id):
    """Extract the setting matrix from the contained scan-dependent crystal
    parameterisations at specified image number."""

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
  def _get_model_data_for_experiment(self, experiment, reflections):
    """helper function to return model data s0, U, B and D for a particular
    experiment. In this scan-varying overload this is trivial because these
    values are already set as arrays in the reflection table"""

    return {'s0':reflections['s0_vector'],
            'U':reflections['u_matrix'],
            'B':reflections['b_matrix'],
            'D':reflections['D_matrix']}

  def _xl_orientation_derivatives(self, isel, parameterisation, reflections):
    """Determine whether dU_dp was precalculated then call the base class
    method"""

    if self._varying_xl_orientations:
      dU_dxlo_p = [reflections["dU_dp{0}".format(i)].select(isel) \
        for i in range(parameterisation.num_free())]
    else:
      dU_dxlo_p = None

    return super(ScanVaryingPredictionParameterisation,
      self)._xl_orientation_derivatives(isel, parameterisation, dU_dxlo_p)

  def _xl_unit_cell_derivatives(self, isel, parameterisation, reflections):
    """Determine whether dB_dp was precalculated then call the base class
    method"""

    if self._varying_xl_unit_cells:
      dB_dxluc_p = [reflections["dB_dp{0}".format(i)].select(isel) \
        for i in range(parameterisation.num_free())]
    else:
      dB_dxluc_p = None

    return super(ScanVaryingPredictionParameterisation,
      self)._xl_unit_cell_derivatives(isel, parameterisation, dB_dxluc_p)

  def calculate_model_state_uncertainties(self, var_cov=None,
                                          obs_image_number=None,
                                          experiment_id=None):
    """Take a variance-covariance matrix of all free parameters (probably
    calculated by a minimisation engine). For each parameterisation in the
    global model, extract the subset of this matrix for the associated block
    of parameters. Pass this on to the relevant model parameterisation to
    calculate its own uncertainty of state.

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
    # identify the crystal parameterisations for this experiment
    xl_op = self._get_xl_orientation_parameterisation(experiment_id)
    xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

    result = {}

    # compose at the requested image number and calculate using the cached
    # varcov matrices. Take the first elt of the list becase the crystal
    # parameterisations are not multi-state
    try:
      xl_op.compose(obs_image_number)
      result['U_cov'] = xl_op.calculate_state_uncertainties(var_cov=None)[0]
    except TypeError:
      pass

    try:
      xl_ucp.compose(obs_image_number)
      result['B_cov'] = xl_ucp.calculate_state_uncertainties(var_cov=None)[0]
    except TypeError:
      pass

    return result

  def set_model_state_uncertainties(self, u_cov_list, b_cov_list,
                                          experiment_id=None):
    """Identify the crystal parameterisations and set the list of covariance
    matrices, if available. They will only be available if the parameterisation
    is a scan-varying type, otherwise they are None"""

    xl_op = self._get_xl_orientation_parameterisation(experiment_id)
    xl_ucp = self._get_xl_unit_cell_parameterisation(experiment_id)

    if u_cov_list:
      try:
        xl_op.set_state_uncertainties(u_cov_list)
      except AttributeError:
        pass

    if b_cov_list:
      try:
        xl_ucp.set_state_uncertainties(b_cov_list)
      except AttributeError:
        pass

    return

class ScanVaryingPredictionParameterisationFast(ScanVaryingPredictionParameterisation):
  """Overloads compose to calculate model states per block rather than per
  reflection"""

  def compose(self, reflections, skip_derivatives=False):
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

      # reset current frame cache for scan-varying parameterisations
      self._current_frame = {}

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

            dmat  = self._get_state_from_parameterisation(dp,
              frame, multi_state_elt=panel_id)
            if dmat is None: dmat = exp.detector[panel_id].get_d_matrix()
            Dmat = exp.detector[panel_id].get_D_matrix()
            reflections['d_matrix'].set_selected(subsel2, dmat)
            reflections['D_matrix'].set_selected(subsel2, Dmat)

            if dp is not None and self._varying_detectors and not skip_derivatives:
              for j, dd in enumerate(dp.get_ds_dp(multi_state_elt=panel_id,
                                                  use_none_as_null=True)):
                if dd is None: continue
                colname = "dd_dp{0}".format(j)
                reflections[colname].set_selected(subsel, dd)

        else: # set states and derivatives for single panel detector

          dmat  = self._get_state_from_parameterisation(dp, frame)
          if dmat is None: dmat = exp.detector[0].get_d_matrix()
          Dmat = exp.detector[0].get_D_matrix()
          reflections['d_matrix'].set_selected(subsel, dmat)
          reflections['D_matrix'].set_selected(subsel, Dmat)

          if dp is not None and self._varying_detectors and not skip_derivatives:
            for j, dd in enumerate(dp.get_ds_dp(use_none_as_null=True)):
              if dd is None: continue
              colname = "dd_dp{0}".format(j)
              reflections[colname].set_selected(subsel, dd)

        # set derivatives of the states for crystal and beam
        if not skip_derivatives:
          if xl_op is not None and self._varying_xl_orientations:
            for j, dU in enumerate(xl_op.get_ds_dp(use_none_as_null=True)):
              if dU is None: continue
              colname = "dU_dp{0}".format(j)
              reflections[colname].set_selected(subsel, dU)
          if xl_ucp is not None and self._varying_xl_unit_cells:
            for j, dB in enumerate(xl_ucp.get_ds_dp(use_none_as_null=True)):
              if dB is None: continue
              colname = "dB_dp{0}".format(j)
              reflections[colname].set_selected(subsel, dB)
          if bp is not None and self._varying_beams:
            for j, ds0 in enumerate(bp.get_ds_dp(use_none_as_null=True)):
              if ds0 is None: continue
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
