#!/usr/bin/env dials.python
#
#
#  Copyright (C) 2014 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Special versions of classes optimised to do detector refinement from multiple
stills in the case where the crystals and beam are fixed."""
from __future__ import division

from dials.algorithms.refinement.refiner import RefinerFactory
from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
  import StillsPredictionParameterisation
from dials.algorithms.refinement.parameterisation.prediction_parameters \
  import SparseGradientVectorMixin
from dials.algorithms.refinement.target_stills import \
  LeastSquaresStillsResidualWithRmsdCutoff
from dials.algorithms.refinement.target import SparseGradientsMixin
from dials.array_family import flex

from dials.algorithms.spot_prediction import ray_intersection

class StillsDetectorRefinerFactory(RefinerFactory):

  @staticmethod
  def config_parameterisation(params, experiments, refman, do_stills):
    """Given a set of parameters, create a parameterisation from a set of
    experimental models.

    Params:
        params The input parameters
        experiments An ExperimentList object

    Returns:
        A tuple of the prediction equation parameterisation and the
        parameter reporter.
    """

    # Shorten parameter paths
    beam_options = params.refinement.parameterisation.beam
    crystal_options = params.refinement.parameterisation.crystal
    detector_options = params.refinement.parameterisation.detector
    sparse = params.refinement.parameterisation.sparse

    # Shorten paths
    import dials.algorithms.refinement.parameterisation as par

    # Parameterise unique Beams
    beam_params = []

    # Parameterise unique Crystals
    xl_ori_params = []
    xl_uc_params = []

    # Parameterise unique Detectors
    det_params = []
    for detector in experiments.detectors():

      exp_ids = experiments.indices(detector)
      # Detector
      if detector_options.panels == "automatic":
        if len(detector) > 1:
          try:
            h = detector.hierarchy()
            det_param = par.DetectorParameterisationHierarchical(detector,
                experiment_ids=exp_ids, level=detector_options.hierarchy_level)
          except AttributeError:
            det_param = par.DetectorParameterisationMultiPanel(detector, beam,
                                                        experiment_ids=exp_ids)
        else:
          det_param = par.DetectorParameterisationSinglePanel(detector,
                                                        experiment_ids=exp_ids)
      elif detector_options.panels == "single":
        det_param = par.DetectorParameterisationSinglePanel(detector,
                                                        experiment_ids=exp_ids)
      elif detector_options.panels == "multiple":
        det_param = par.DetectorParameterisationMultiPanel(detector, beam,
                                                        experiment_ids=exp_ids)
      elif detector_options.panels == "hierarchical":
        det_param = par.DetectorParameterisationHierarchical(detector, beam,
                experiment_ids=exp_ids, level=detector_options.hierarchy_level)
      else: # can only get here if refinement.phil is broken
        raise RuntimeError("detector_options.panels value not recognised")

      if detector_options.fix:
        if detector_options.fix == "all":
          det_param.set_fixed([True] * det_param.num_total())
        elif detector_options.fix == "position":
          to_fix = [e.param_type.startswith('length') \
                    for e in det_param.get_params(only_free = False)]
          det_param.set_fixed(to_fix)
        elif detector_options.fix == "orientation":
          to_fix = [e.param_type.startswith('angle') \
                    for e in det_param.get_params(only_free = False)]
          det_param.set_fixed(to_fix)
        else: # can only get here if refinement.phil is broken
          raise RuntimeError("detector_options.fix value not recognised")

      if detector_options.fix_list:
        to_fix = [True if i in detector_options.fix_list else False \
                  for i in range(det_param.num_total())]
        det_param.set_fixed(to_fix)

      det_params.append(det_param)

    # Now we have the final list of model parameterisations, build a restraints
    # parameterisation (if requested). Only unit cell restraints are supported
    # at the moment.
    if any([crystal_options.unit_cell.restraints.tie_to_target,
            crystal_options.unit_cell.restraints.tie_to_group]):
      restraints_param = cls.config_restraints(params, det_params, beam_params,
        xl_ori_params, xl_uc_params)
    else:
      restraints_param = None

    # Prediction equation parameterisation
    if do_stills: # doing stills
      if sparse:
        spp = StillsDetectorPredictionParameterisationSparse
      else:
        spp = StillsDetectorPredictionParameterisation
      pred_param = spp(experiments, det_params, beam_params,
                       xl_ori_params, xl_uc_params)

    else: # doing scans
      raise NotImplementedError("currently only for stills")
      #if crystal_options.scan_varying:
      #  if crystal_options.UB_model_per == "reflection":
      #    #from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
      #    #  import ScanVaryingPredictionParameterisation as PredParam
      #    raise NotImplementedError("currently only for stills")
      #  elif crystal_options.UB_model_per == "image":
      #    #from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
      #    #  import ScanVaryingPredictionParameterisationFast as PredParam
      #    raise NotImplementedError("currently only for stills")
      #  else:
      #    raise RuntimeError("UB_model_per=" + crystal_options.scan_varying +
      #                       " is not a recognised option")
      #  pred_param = PredParam(
      #        experiments,
      #        det_params, beam_params, xl_ori_params, xl_uc_params)
      #else:
      #  if sparse:
      #    #from dials.algorithms.refinement.parameterisation.prediction_parameters \
      #    #  import XYPhiPredictionParameterisationSparse as PredParam
      #    raise NotImplementedError("currently only for stills")
      #  else:
      #    #from dials.algorithms.refinement.parameterisation.prediction_parameters \
      #    #  import XYPhiPredictionParameterisation as PredParam
      #    raise NotImplementedError("currently only for stills")
      #  pred_param = PredParam(
      #      experiments,
      #      det_params, beam_params, xl_ori_params, xl_uc_params)

    # Parameter reporting
    param_reporter = par.ParameterReporter(det_params, beam_params,
                                           xl_ori_params, xl_uc_params)

    return pred_param, param_reporter, restraints_param

  @staticmethod
  def config_target(params, experiments, refman, do_stills):
    """Given a set of parameters, configure a factory to build a
    target function

    Params:
        params The input parameters

    Returns:
        The target factory instance
    """

    # Shorten parameter paths
    options = params.refinement.target
    sparse = params.refinement.parameterisation.sparse

    if options.rmsd_cutoff == "fraction_of_bin_size":
      absolute_cutoffs = None
    elif options.rmsd_cutoff == "absolute":
      absolute_cutoffs = options.absolute_cutoffs
    else:
      raise RuntimeError("Target function rmsd_cutoff option" +
          options.rmsd_cutoff + " not recognised")

    # all experiments have the same (or no) goniometer
    goniometer = experiments[0].goniometer
    for e in experiments: assert e.goniometer is goniometer

    # build managed reflection predictors
    from dials.algorithms.refinement.prediction import ExperimentsPredictor
    ref_predictor = ExperimentsPredictor(experiments, do_stills)

    # Determine whether the target is in X, Y, Phi space or just X, Y.
    if do_stills:
      if sparse:
        targ = LeastSquaresStillsDetectorSparse
      else:
        targ = LeastSquaresStillsDetector
    else:
      raise NotImplementedError("currently only for stills")
      #if sparse:
      #  raise NotImplementedError("currently only for stills")
      #  #from dials.algorithms.refinement.target \
      #  #  import LeastSquaresPositionalResidualWithRmsdCutoffSparse as targ
      #else:
      #  raise NotImplementedError("currently only for stills")
      #  #from dials.algorithms.refinement.target \
      #  #  import LeastSquaresPositionalResidualWithRmsdCutoff as targ

    # Here we pass in None for prediction_parameterisation and
    # restraints_parameterisation, as these will be linked to the object later
    target = targ(experiments=experiments,
                  reflection_predictor=ref_predictor,
                  ref_man=refman,
                  prediction_parameterisation=None,
                  restraints_parameterisation=None,
                  frac_binsize_cutoff=options.bin_size_fraction,
                  absolute_cutoffs=absolute_cutoffs,
                  gradient_calculation_blocksize=options.gradient_calculation_blocksize)

    return target


class StillsDetectorPredictionParameterisation(StillsPredictionParameterisation):

  def get_gradients(self, reflections, callback=None):
    """
    Calculate gradients of the prediction formula with respect to each
    of the parameters of the detector, for all of the reflections.

    """

    ### Calculate various quantities of interest for the reflections

    # Set up arrays of values for each reflection
    n = len(reflections)
    D = flex.mat3_double(n)
    #s0 = flex.vec3_double(n)
    #U = flex.mat3_double(n)
    #B = flex.mat3_double(n)
    #axis = flex.vec3_double(n)

    for iexp, exp in enumerate(self._experiments):

      sel = reflections['id'] == iexp
      isel = sel.iselection()

      # D matrix array
      panels = reflections['panel'].select(isel)
      for ipanel, D_mat in enumerate([p.get_D_matrix() for p in exp.detector]):
        subsel = isel.select(panels == ipanel)
        D.set_selected(subsel, D_mat)

      # s0 array
      #s0.set_selected(isel, exp.beam.get_s0())

      # U and B arrays
      #exp_U, exp_B = self._get_U_B_for_experiment(exp.crystal, reflections, isel)
      #U.set_selected(isel, exp_U)
      #B.set_selected(isel, exp_B)

      # axis array
      #if exp.goniometer:
      #  axis.set_selected(isel, exp.goniometer.get_rotation_axis())
    return self._get_gradients_core(reflections, D, callback)

  def _get_gradients_core(self, reflections, D, callback=None):
    """Calculate gradients of the prediction formula with respect to
    each of the parameters of the contained models, for reflection h
    with scattering vector s that intersects panel panel_id. That is,
    calculate dX/dp, dY/dp and dDeltaPsi/dp. Ignore axis because these
    are stills"""

    # pv is the 'projection vector' for the ray along s1.
    self._D = D
    self._s1 = reflections['s1']
    self._pv = D * self._s1

    # also need quantities derived from pv, precalculated for efficiency
    u, v, w = self._pv.parts()
    self._w_inv = 1/w
    self._u_w_inv = u * self._w_inv
    self._v_w_inv = v * self._w_inv

    self._DeltaPsi = reflections['delpsical.rad']

    # q is the reciprocal lattice vector, in the lab frame
    #self._UB = U * B
    #self._U = U
    #self._B = B
    #self._h = reflections['miller_index'].as_vec3_double()
    #self._q = (self._UB * self._h)
    #self._q_scalar = self._q.norms()
    #self._qq = self._q_scalar * self._q_scalar

    # r is the reciprocal lattice vector rotated to the Ewald sphere
    #self._s0 = s0
    #self._r = self._s1 - self._s0

    # we also need the unit directions q0 and s0u
    #self._q0 = self._q.each_normalize()
    #self._s0u = self._s0.each_normalize()

    # e1 is the unit vector about which DeltaPsi rotation is defined
    #self._e1 = self._q0.cross(self._s0u).each_normalize()

    # q1 completes an orthonormal set with q0 and e1
    #self._q1 = self._q0.cross(self._e1).each_normalize()

    # we want the wavelength
    #self._wavelength = 1. / self._s0.norms()

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
    ### to derivatives d[pv]/dp and d[DeltaPsi]/dp

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
      for panel_id, _ in enumerate(detector):

        # get the right subset of array indices to set for this panel
        sub_isel = isel.select(panel == panel_id)
        if len(sub_isel) == 0:
          # if no reflections intersect this panel, skip calculation
          continue
        sub_pv = self._pv.select(sub_isel)
        sub_D = self._D.select(sub_isel)
        dpv_ddet_p = self._detector_derivatives(dp, sub_pv, sub_D, panel_id)

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
          results[iparam][self._grad_names[0]].set_selected(sub_isel, dX)
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

    return results

class StillsDetectorPredictionParameterisationSparse(SparseGradientVectorMixin,
  StillsDetectorPredictionParameterisation):
  pass


class LeastSquaresStillsDetector(LeastSquaresStillsResidualWithRmsdCutoff):

  _first_predict = True
  def predict(self):
    """perform reflection prediction and update the reflection manager"""

    if self._first_predict:
      self._first_predict = False
      super(LeastSquaresStillsDetector, self).predict()

      # HACK TO PUT IN A PHI COLUMN, WHICH RAY_INTERSECTION EXPECTS
      reflections = self._reflection_manager.get_obs()
      reflections['phi'] = flex.double(len(reflections), 0)
      return
    else:

      # update the reflection_predictor with the scan-independent part of the
      # current geometry
      #self._reflection_predictor.update()

      # reset the 'use' flag for all observations
      #self._reflection_manager.reset_accepted_reflections()

      # do prediction (updates reflection table in situ).
      reflections = self._reflection_manager.get_obs()
      #self._reflection_predictor.predict(reflections)
      # FIXME HACK TO GET THE DETECTOR FROM THE FIRST EXPERIMENT
      detector = self._reflection_predictor._experiments[0].detector
      success = ray_intersection(detector, reflections, reflections['panel'])
      assert success.all_eq(True)

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
      #mask = reflections.get_flags(reflections.flags.predicted)
      #reflections.set_flags(mask, reflections.flags.used_in_refinement)

      # collect the matches
      self.update_matches(force=True)

    return

class LeastSquaresStillsDetectorSparse(SparseGradientsMixin,
  LeastSquaresStillsDetector):
  pass
