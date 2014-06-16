#!/usr/bin/env python
#
# dials.algorithms.refinement.refiner.py
#
#  Copyright (C) 2013 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Authors: James Parkhurst, David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Refiner is the refinement module public interface. RefinerFactory is
what should usually be used to construct a Refiner."""

from __future__ import division
from dials.algorithms.refinement.refinement_helpers import print_model_geometry
from dials.model.experiment.experiment_list import ExperimentList, Experiment
from dials.array_family import flex # import dependency

class RefinerFactory(object):
  """Factory class to create refiners"""

  @classmethod
  def from_parameters_data_experiments(cls,
                                       params,
                                       reflections,
                                       experiments,
                                       verbosity=None):

    #TODO Checks on the input
    #E.g. does every experiment contain at least one overlapping model with at
    #least one other experiment? Are all the experiments either rotation series
    #or stills (the combination of both not yet supported)?

    # if no verbosity override is given, take from the parameters
    if not verbosity:
      verbosity = params.refinement.verbosity

    # copy the experiments
    import copy
    experiments = copy.deepcopy(experiments)

    # copy the reflections
    reflections = copy.deepcopy(reflections)

    return cls._build_components(params,
                                 reflections,
                                 experiments,
                                 crystal_ids=range(len(experiments)),
                                 verbosity=verbosity)

  @classmethod
  def from_parameters_data_models(cls,
                                  params,
                                  reflections,
                                  sweep=None,
                                  beam=None,
                                  goniometer=None,
                                  detector=None,
                                  scan=None,
                                  image_width_rad=None,
                                  sweep_range_rad=None,
                                  crystal=None,
                                  crystals=None,
                                  crystal_ids=None,
                                  verbosity=None):
    """Given a set of parameters, reflections and experimental models for a
    single Experiment, construct a refiner.

    Mandatory arguments:
      params - The input parameters as a phil scope_extract object
      reflections - Input ReflectionList data

    Argument alternatives:
      sweep - Experimental models as a dxtbx sweep
        or
      /beam - A dxtbx Beam object
      \detector - A dxtbx Detector object
      crystal - A cctbx crystal_model object
        or
      /crystals - A list of cctbx crystal_model objects
      \crystal_ids - A list of integer crystal ids to match to crystals

    Optional arguments:
      goniometer - A dxtbx Goniometer model
      scan - A dxtbx Scan model
        or
      /image_width_rad - The 'width' of one image in radians
      \sweep_range_rad - Pair of rotation scan extremes in radians
      verbosity - An integer verbosity level

    Returns:
      The Refiner instance.

    Notes
    -----

    The interface is intended to be flexible by allowing input to be passed in
    various forms. Some are incompatible, e.g. passing a sweep alongside any
    other dxtbx model is disallowed. Either a crystal or a list of crystals must
    be provided (but not both). A list of crystals must be associated with a
    list of crystal id integers to associate with the crystals.

    The optional arguments determine the behaviour of the refiner, alongside the
    phil parameters. Of particular interest, the presence of a goniometer model
    decides whether the refinement is of a rotation scan, with a target
    expressed in (X, Y, Phi) space, or of a still image, with an (X, Y) target.
    For a rotation scan, a Scan object is optional. It is only used for some
    information, which may be provided instead by using the image_width_rad and
    sweep_range_rad arguments. It is even possible to specify none of
    image_width_rad, sweep_range_rad or scan and still do (X, Y, Phi)
    refinement, but in such a case the phil params should not specify
    reflections.reflections_per_degree and must specify
    target.rmsd_cutoff=absolute. This is checked to avoid illegal construction
    of the target and reflection manager objects.

    The steps performed by this factory function are:
    * check for consistent input
    * copy the input models and reflections
    * set beam vectors in the working copy of reflections if not already set
    * create parameterisations of the models that were provided, depending on
      phil preferences for choices such as fixing certain parameters
    * create a reflection manager, depending on phil preferences for decisions
      regarding inclusion and random sampling criteria
    * create a target function using phil parameters to determine the 'target
      achieved' criterion, and the presence of a goniometer to determine whether
      to include or ignore residuals in phi
    * create a refinery (minimisation engine) using phil parameters to select a
      particular algorithm, determine a maximum number of steps and to control
      additional logging features
    * package all the above together in a 'refiner' object, which provides the
      interface for running refinement
    * return that refiner

    """

    # if no verbosity override is given, take from the parameters
    if not verbosity:
      verbosity = params.refinement.verbosity

    # checks on the input
    if sweep:
      assert [beam, goniometer, detector, scan].count(None) == 4
      beam = sweep.get_beam()
      detector = sweep.get_detector()
      goniometer = sweep.get_goniometer()
      scan = sweep.get_scan()

    # if either of these provided, both must be
    if image_width_rad or sweep_range_rad:
      assert image_width_rad and sweep_range_rad

    # only one of image_width_rad or scan should be provided (or neither)
    assert [image_width_rad, scan].count(None) >= 1
    #if scan:
      #image_width_rad = scan.get_oscillation(deg=False)[1]
      #sweep_range_rad = scan.get_oscillation_range(deg=False)
    if image_width_rad:
      # create a dummy scan
      from dxtbx.model.scan import scan_factory
      sf = scan_factory()
      sweep_width = abs(sweep_range_rad[1] - sweep_range_rad[0])
      nimages = int(sweep_width / image_width_rad)
      scan = sf.make_scan(image_range= (1, nimages),
                          exposure_times = 0.,
                          oscillation = (0, image_width_rad),
                          epochs = range(nimages),
                          deg = False)

    if goniometer and not image_width_rad and not scan:
      # This is not to be supported any more. Now we state we must have enough
      # information to make a scan
      raise RuntimeError("Goniometer provided, but not enough information about a scan!")
      # if there was neither a scan nor the image width provided,
      # the target rmsd must be provided in absolute terms
      #assert params.refinement.target.rmsd_cutoff == "absolute"
      # also there is no sweep range, so the reflection manager
      # cannot sample by number of reflections per degree
      #assert params.refinement.reflections.reflections_per_degree is None

    # do we have the essential models?
    assert [beam, detector].count(None) == 0
    assert [crystal, crystals].count(None) == 1
    if crystals: assert len(crystals) == len(crystal_ids)

    # copy the models
    from dxtbx.model import Beam
    import copy
    # use copy constructor
    beam = Beam(beam)
    detector = copy.deepcopy(detector)
    if crystal:
      crystals = [copy.deepcopy(crystal)]
      crystal_ids = [0]
    if crystals: crystals = copy.deepcopy(crystals)

    # build an experiment list with the copied models
    experiments = ExperimentList()
    for crystal in crystals:
      experiments.append(Experiment(
        beam=beam, detector=detector, goniometer=goniometer,
        scan=scan, crystal=crystal, imageset=None))

    # copy the reflections
    reflections = copy.deepcopy(reflections)

    # Build components and return
    return cls._build_components(params,
                                 reflections,
                                 experiments,
                                 crystal_ids,
                                 verbosity)

  @classmethod
  def _build_components(cls, params, reflections, experiments, crystal_ids,
                        verbosity):
    """low level build"""

    # check that the beam vectors are stored: if not, compute them
    refs_wo_s1_sel = (reflections['s1'].norms() < 1.e-6)
    nrefs_wo_s1 = refs_wo_s1_sel.count(True)
    if nrefs_wo_s1 > 0 and verbosity > 1:
      print "Setting scattering vectors for", nrefs_wo_s1, "reflections"
    for i_expt, expt in enumerate(experiments):
      detector = expt.detector
      beam = expt.beam
      expt_sel = reflections['id'] == i_expt
      for i_panel, panel in enumerate(detector):
        panel_sel = reflections['panel'] == i_panel
        isel = (expt_sel & panel_sel & refs_wo_s1_sel).iselection()
        spots = reflections.select(isel)
        x, y, rot_angle = spots['xyzobs.mm.value'].parts()
        s1 = panel.get_lab_coord(flex.vec2_double(x,y))
        s1 = s1/s1.norms() * (1/beam.get_wavelength())
        reflections['s1'].set_selected(isel, s1)

    # unset the refinement flags (creates flags field if needed)
    from dials.array_family.flex import reflection_table
    reflections.unset_flags(flex.size_t_range(len(reflections)),
        reflection_table.flags.used_in_refinement)

    # create parameterisations
    pred_param, param_reporter = \
            cls.config_parameterisation(params, experiments)

    if verbosity > 1:
      print "Prediction equation parameterisation built\n"
      print "Parameter order : name mapping"
      for i, e in enumerate(pred_param.get_param_names()):
        print "Parameter %03d : " % i + e
      print

    if verbosity > 1:
      print "Building reflection manager"
      print ("Input reflection list size = %d observations"
             % len(reflections))

    # create reflection manager
    refman = cls.config_refman(params, reflections, experiments, verbosity)

    if verbosity > 1:
      print ("Number of observations that pass initial inclusion criteria = %d"
             % refman.get_accepted_refs_size())
      sample_size = refman.get_sample_size()
      if sample_size:
        print ("Working set size = %d observations" % sample_size)
      print "Reflection manager built\n"

    if verbosity > 1: print "Building target function"

    # create target function
    target = cls.config_target(params, experiments, refman, pred_param)

    if verbosity > 1: print "Target function built\n"

    if verbosity > 1: print "Building refinement engine"

    # create refinery
    refinery = cls.config_refinery(params, target, pred_param, verbosity)

    if verbosity > 1: print "Refinement engine built\n"

    # build refiner interface and return
    return Refiner(reflections, experiments, crystal_ids,
                    pred_param, param_reporter, refman, target, refinery,
                    verbosity=verbosity)

  @staticmethod
  def config_parameterisation(params, experiments):
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

    # Currently a refinement job can only have one parameterisation of the
    # prediction equation. This can either be of the XYDelPsi (stills) type, the
    # XYPhi (scans) type or the scan-varying XYPhi type with a varying crystal
    # model

    # Parameterise unique Beams
    beam_params_stills = []
    beam_params_scans = []
    for beam in experiments.beams():
      # The Beam is parameterised with reference to a goniometer axis (or None).
      # Therefore this Beam must always be found alongside the same Goniometer
      # (otherwise it should be a different Beam as it will require a different
      # parameterisation).
      exp_ids = experiments.indices(beam)
      assoc_gonios = [experiments[i].goniometer for i in exp_ids]
      goniometer = assoc_gonios[0]
      assert all(g is goniometer for g in assoc_gonios)

      # Parameterise, passing the goniometer (but accepts None)
      beam_param = par.BeamParameterisationOrientation(beam, goniometer,
                                                       experiment_ids=exp_ids)
      if beam_options.fix:
        if beam_options.fix == "all":
          beam_param.set_fixed([True, True])
        elif beam_options.fix == "in_spindle_plane":
          beam_param.set_fixed([True, False])
        elif beam_options.fix == "out_spindle_plane":
          beam_param.set_fixed([False, True])
        else: # can only get here if refinement.phil is broken
          raise RuntimeError("beam_options.fix value not recognised")

      if beam_options.fix_list:
        to_fix = [True if i in beam_options.fix_list else False \
                  for i in range(beam_param.num_total())]
        beam_param.set_fixed(to_fix)

      if goniometer is None:
        beam_params_stills.append(beam_param)
      else:
        beam_params_scans.append(beam_param)

    # Parameterise unique Crystals
    xl_ori_params_stills = []
    xl_ori_params_scans = []
    xl_uc_params_stills = []
    xl_uc_params_scans = []
    for crystal in experiments.crystals():
      # This crystal can only ever appear either in scans or in stills
      # (otherwise it requires a different crystal model)
      exp_ids = experiments.indices(crystal)
      assoc_models = [(experiments[i].goniometer, experiments[i].scan) \
                      for i in exp_ids]
      goniometer, scan = assoc_models[0]
      if goniometer is None:
        assert all(g is None and s is None for (g, s) in assoc_models)

      if crystal_options.scan_varying:
        # If a crystal is scan-varying, then it must always be found alongside
        # the same Scan and Goniometer in any Experiments in which it appears
        assert [goniometer, scan].count(None) == 0
        assert all(g is goniometer and s is scan for (g, s) in assoc_models)

        if crystal_options.num_intervals == "fixed_width":
          sweep_range_deg = scan.get_oscillation_range(deg=True)
          deg_per_interval = crystal_options.interval_width_degrees
          n_intervals = int(
            abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval)
        else:
          n_intervals = crystal_options.absolute_num_intervals

        xl_ori_param = par.ScanVaryingCrystalOrientationParameterisation(
                                            crystal,
                                            scan.get_array_range(),
                                            n_intervals,
                                            experiment_ids=exp_ids)
        xl_uc_param = par.ScanVaryingCrystalUnitCellParameterisation(
                                            crystal,
                                            scan.get_array_range(),
                                            n_intervals,
                                            experiment_ids=exp_ids)
      else:
        xl_ori_param = par.CrystalOrientationParameterisation(crystal,
                                                        experiment_ids=exp_ids)
        xl_uc_param = par.CrystalUnitCellParameterisation(crystal,
                                                        experiment_ids=exp_ids)

      if crystal_options.fix:
        if crystal_options.fix == "all":
          xl_ori_param.set_fixed([True] * xl_ori_param.num_total())
          xl_uc_param.set_fixed([True] * xl_uc_param.num_total())
        elif crystal_options.fix == "cell":
          xl_uc_param.set_fixed([True] * xl_uc_param.num_total())
        elif crystal_options.fix == "orientation":
          xl_ori_param.set_fixed([True] * xl_ori_param.num_total())
        else: # can only get here if refinement.phil is broken
          raise RuntimeError("crystal_options.fix value not recognised")

      if crystal_options.cell_fix_list:
        to_fix = [True if i in crystal_options.cell_fix_list else False \
                  for i in range(xl_uc_param.num_total())]
        xl_uc_param.set_fixed(to_fix)

      if crystal_options.orientation_fix_list:
        to_fix = [True if i in crystal_options.orientation_fix_list else False \
                  for i in range(xl_ori_param.num_total())]
        xl_ori_param.set_fixed(to_fix)

      if goniometer is None:
        xl_ori_params_stills.append(xl_ori_param)
        xl_uc_params_stills.append(xl_uc_param)
      else:
        xl_ori_params_scans.append(xl_ori_param)
        xl_uc_params_scans.append(xl_uc_param)

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

    # Determine whether this is a parameterisation for scans or for stills
    if beam_params_stills == []:
      # we expect scans
      assert xl_ori_params_stills == []
      assert xl_uc_params_stills == []
      assert beam_params_scans != []
      assert xl_ori_params_scans != []
      assert xl_uc_params_scans != []
      beam_params = beam_params_scans
      xl_ori_params = xl_ori_params_scans
      xl_uc_params = xl_uc_params_scans
      param_type = "scans"
    else:
      # we expect stills
      assert beam_params_scans == []
      assert xl_ori_params_scans == []
      assert xl_uc_params_scans == []
      assert xl_ori_params_stills != []
      assert xl_uc_params_stills != []
      beam_params = beam_params_stills
      xl_ori_params = xl_ori_params_stills
      xl_uc_params = xl_uc_params_stills
      param_type = "stills"

    # Prediction equation parameterisation
    if param_type is "scans":
      if crystal_options.scan_varying:
        if crystal_options.UB_model_per == "reflection":
          from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
            import VaryingCrystalPredictionParameterisation as PredParam
        elif crystal_options.UB_model_per == "image":
          from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
            import VaryingCrystalPredictionParameterisationFast as PredParam
        else:
          raise RuntimeError("UB_model_per=" + crystal_options.scan_varying +
                             " is not a recognised option")
        pred_param = PredParam(
              experiments,
              det_params, beam_params, xl_ori_params, xl_uc_params)
      else:
        # FIXME temporary user choice for refactored version
        if sparse:
          from dials.algorithms.refinement.parameterisation.prediction_parameters \
            import XYPhiPredictionParameterisationDebug as PredParam
        else:
          from dials.algorithms.refinement.parameterisation.prediction_parameters \
            import XYPhiPredictionParameterisation as PredParam
        pred_param = PredParam(
            experiments,
            det_params, beam_params_scans, xl_ori_params_scans, xl_uc_params_scans)
    else:
      assert param_type is "stills"
      # FIXME temporary user choice for sparse matrix version
      # currently only using the sparse version for stills.
      if sparse:
        from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import StillsPredictionParameterisationSparse as StillsPredictionParameterisation
      else:
        from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import StillsPredictionParameterisation
      pred_param = StillsPredictionParameterisation(
          experiments,
          det_params, beam_params, xl_ori_params, xl_uc_params)

    # Parameter reporting
    param_reporter = par.ParameterReporter(det_params, beam_params,
                                           xl_ori_params, xl_uc_params)

    return pred_param, param_reporter

  @staticmethod
  def config_refinery(params, target, pred_param, verbosity):
    """Given a set of parameters, a target class, and a prediction
    parameterisation class, build a refinery

    Params:
        params The input parameters

    Returns:
        The refinery instance
    """

    # Shorten parameter path
    options = params.refinement.refinery

    import dials.algorithms.refinement.engine as engine

    if options.engine == "SimpleLBFGS":
      from engine import SimpleLBFGS as refinery
    elif options.engine == "LBFGScurvs":
      from engine import LBFGScurvs as refinery
    elif options.engine == "GaussNewton":
      from engine import GaussNewtonIterations as refinery
    elif options.engine == "LevMar":
      from engine import LevenbergMarquardtIterations as refinery
    else:
      raise RuntimeError("Refinement engine " + options.engine +
                         " not recognised")

    if verbosity > 1: print "Selected refinement engine type:", options.engine

    return refinery(target = target,
            prediction_parameterisation = pred_param,
            log = options.log,
            verbosity = verbosity,
            track_step = options.track_step,
            track_gradient = options.track_gradient,
            track_parameter_correlation = options.track_parameter_correlation,
            max_iterations = options.max_iterations)

  @staticmethod
  def config_refman(params, reflections, experiments, verbosity):
    """Given a set of parameters and models, build a reflection manager

    Params:
        params The input parameters

    Returns:
        The reflection manager instance
    """

    # Shorten parameter path
    options = params.refinement.reflections
    options.random_seed
    if options.use_all_reflections:
      nref_per_degree = None
    else:
      nref_per_degree = options.reflections_per_degree

    # While a random subset of reflections is used, continue to
    # set random.seed to get consistent behaviour
    if options.random_seed:
      import random
      random.seed(options.random_seed)
      flex.set_random_seed(options.random_seed)
      if verbosity > 1:
        print "Random seed set to %d\n" % options.random_seed

    if all(e.goniometer is not None for e in experiments):
      from dials.algorithms.refinement.reflection_manager import ReflectionManager as refman
    elif all(e.goniometer is None for e in experiments):
      from dials.algorithms.refinement.reflection_manager import \
          StillsReflectionManager as refman

    else:
      raise NotImplementedError("ExperimentList contains a mixture of "
        "experiments with goniometers and those without. This is not currently "
        "supported.")

    # do outlier rejection?
    if options.do_outlier_rejection:
      iqr_multiplier=options.iqr_multiplier
    else:
      iqr_multiplier=None

    return refman(reflections=reflections,
            experiments=experiments,
            nref_per_degree=nref_per_degree,
            min_num_obs=options.minimum_number_of_reflections,
            max_sample_size = options.maximum_sample_size,
            min_sample_size = options.minimum_sample_size,
            close_to_spindle_cutoff=options.close_to_spindle_cutoff,
            iqr_multiplier=iqr_multiplier,
            verbosity=verbosity)

  @staticmethod
  def config_target(params, experiments, refman, pred_param):
    """Given a set of parameters, configure a factory to build a
    target function

    Params:
        params The input parameters

    Returns:
        The target factory instance
    """

    # Shorten parameter path
    options = params.refinement.target

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
    ref_predictor = ExperimentsPredictor(experiments)

    # Determine whether the target is in X, Y, Phi space or just X, Y.
    if all(e.goniometer is not None for e in experiments):
      from dials.algorithms.refinement.target \
        import LeastSquaresPositionalResidualWithRmsdCutoff as targ
    elif all(e.goniometer is None for e in experiments):
      from dials.algorithms.refinement.target_stills \
        import LeastSquaresStillsResidualWithRmsdCutoff as targ
    else:
      raise NotImplementedError("ExperimentList contains a mixture of "
        "experiments with goniometers and those without. This is not currently "
        "supported.")

    target = targ(experiments, ref_predictor, refman, pred_param,
                    options.bin_size_fraction, absolute_cutoffs,
                    options.jacobian_max_nref)

    return target

class Refiner(object):
  """Public interface for performing DIALS refinement.

  Public methods:
    run
    rmsds
    get_experiments
    get_beam
    get_crystal
    get_crystals
    get_detector
    get_goniometer
    get_scan
    get_reflections
    get_matches
    get_param_reporter
    parameter_correlation_plot
    selection_used_for_refinement
    predict_for_reflection_table

  Notes:
    * The return value of run is a recorded history of the refinement
    * The experiments accessor provides a copy of the experiments used by
      refinement
    * The model accessors provide copies of those models that might be modified
      by refinement (beam, crystal and detector) TO BE DEPRECATED
    * get_matches exposes the function of the same name from the privately
      stored reflection manager
    * The return value of selection_used_for_refinement is a flex.bool
    * predict_reflections does so for all observable reflections, not just the
      set passed into the refiner

    """

  def __init__(self, reflections, experiments, crystal_ids,
               pred_param, param_reporter, refman, target, refinery,
               verbosity=0):
    """
    Mandatory arguments:
      reflections - Input ReflectionList data
      beam - A dxtbx Beam object
      crystals - A list of cctbx crystal_model objects
      detector - A dxtbx Detector object
      pred_param - An object derived from the PredictionParameterisation class
      param_reporter -A ParameterReporter object
      refman - A ReflectionManager object
      target - An object derived from the Target class
      refinery - An object derived from the Refinery class

    Optional arguments:
      goniometer - A dxtbx Goniometer object
      scan - A dxtbx Scan object
      verbosity - An integer verbosity level

    """

    # FIXME only models from the first Experiment are kept here
    detector = experiments[0].detector
    beam = experiments[0].beam
    crystals = experiments.crystals()
    goniometer = experiments[0].goniometer
    scan = experiments[0].scan

    # keep the data and models public for access after refinement
    self._experiments = experiments
    self._reflections = reflections
    self._beam = beam
    self._crystals = crystals
    self._crystal_ids = crystal_ids
    # only keep crystal if there is indeed only one
    self._crystal = crystals[0] if len(crystals) == 1 else None
    self._detector = detector

    # these could be None (for stills/XFEL)
    self._goniometer = goniometer
    self._scan = scan

    # refinement module main objects
    self._pred_param = pred_param
    self._refman = refman
    self._target = target
    self._refinery = refinery

    # parameter reporter
    self._param_report = param_reporter

    self._verbosity = verbosity

    return

  def get_experiments(self):

    from copy import deepcopy
    return deepcopy(self._experiments)

  def get_beam(self):
    """Return a copy of the beam model"""

    from dxtbx.model import Beam
    return Beam(self._beam)

  def get_crystal(self):
    """Return a copy of the crystal model, if present"""

    from copy import deepcopy
    if self._crystal:
      return deepcopy(self._crystal)
    else:
      return None

  def get_crystals(self):
    """Return a copy of the crystal models (always present)"""

    from copy import deepcopy
    return deepcopy(self._crystals)

  def get_detector(self):
    """Return a copy of the detector model"""

    from copy import deepcopy
    from dxtbx.model import Detector
    return deepcopy(self._detector)

  def get_goniometer(self):
    """Return the goniometer model, if present"""

    # this is unmodified by refinement, so safe to return directly without copy
    if self._goniometer:
      return self._goniometer
    else:
      return None

  def get_scan(self):
    """Return the scan model, if present"""

    # this is unmodified by refinement, so safe to return directly without copy
    if self._scan:
      return self._scan
    else:
      return None

  def get_reflections(self):
    """Return the input reflections"""

    # FIXME Consider: Does a Refiner really need to keep a copy of the input
    # reflections? (indexing code seems to use it, but is it necessary?)

    # The length and order of this is unmodified by refinement - only beam
    # vectors may have been set, if they were not already. It is deemed safe
    # to return this without a copy (NB the input will usually already have been
    # copied by the RefinerFactory)
    return self._reflections

  def rmsds(self):
    """Return rmsds of the current model"""

    self._refinery.prepare_for_step()

    return self._target.rmsds()

  def get_matches(self):
    """Delegated to the reflection manager"""

    # FIXME Consider: Does this information really need to be exposed by the
    # public API (indexing code seems to use it, but is it necessary?)
    return self._refman.get_matches()

  def get_param_reporter(self):
    """Get the ParameterReport object linked to this Refiner"""

    return self._param_report

  def parameter_correlation_plot(self, step):
    """Create a correlation matrix plot between columns of the Jacobian at
    the specified refinement step. Inspired by R's corrplot and
    https://github.com/louridas/corrplot/blob/master/corrplot.py"""

    corrmat = self._refinery.get_correlation_matrix_for_step(step)
    if corrmat is None: return None

    from math import pi, sqrt
    try:
      import matplotlib.pyplot as plt
      import matplotlib.cm as cm
    except ImportError as e:
      print "matplotlib modules not available", e
      return None

    labels = self._pred_param.get_param_names()
    plt.figure(1)
    ax = plt.subplot(1, 1, 1, aspect='equal')
    width, height = corrmat.all()
    num_cols, num_rows = width, height
    poscm = cm.get_cmap('Blues')
    negcm = cm.get_cmap('Reds')
    for x in xrange(width):
      for y in xrange(height):
        d = corrmat[x, y]
        rotate = -45 if d > 0 else +45
        clrmap = poscm if d >= 0 else negcm
        d_abs = abs(d)
        circ = plt.Circle((x, y),radius=sqrt(d_abs/pi))
        circ.set_edgecolor('white')
        circ.set_facecolor(clrmap(d_abs))
        ax.add_artist(circ)
    ax.set_xlim(-1, num_cols)
    ax.set_ylim(-1, num_rows)

    ax.xaxis.tick_top()
    xtickslocs = range(len(labels))
    ax.set_xticks(xtickslocs)
    ax.set_xticklabels(labels, rotation=30, fontsize='small', ha='left')

    ax.invert_yaxis()
    ytickslocs = range(len(labels))
    ax.set_yticks(ytickslocs)
    ax.set_yticklabels(labels, fontsize='small')

    return plt

  def run(self):
    """Run refinement"""

    ####################################
    # Do refinement and return history #
    ####################################

    if self._verbosity > 1:
      print ""
      print "Experimental models before refinement"
      print "-------------------------------------"
      print self._beam
      print self._detector
      if self._goniometer: print self._goniometer
      if self._scan: print self._scan
      for i, x in zip(self._crystal_ids, self._crystals): print i, x

    if self._verbosity > 0:
      print ""
      print "Running refinement"
      print "------------------"
    self._refinery.run()

    # write scan varying setting matrices back to crystal models
    #FIXME tidy up
    from dials.algorithms.refinement.parameterisation import \
      VaryingCrystalPredictionParameterisation
    from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
      import VaryingCrystalPredictionParameterisationFast
    if isinstance(self._pred_param, VaryingCrystalPredictionParameterisation) or \
       isinstance(self._pred_param, VaryingCrystalPredictionParameterisationFast):
      for iexp, exp in enumerate(self._experiments):
        ar_range = exp.scan.get_array_range()
        A_list = [self._pred_param.get_UB(t, iexp) for t in range(ar_range[0],
                                                            ar_range[1]+1)]
        exp.crystal.set_A_at_scan_points(A_list)

    if self._verbosity > 1:
      print
      print "Experimental models after refinement"
      print "------------------------------------"
      print self._beam
      print self._detector
      if self._goniometer: print self._goniometer
      if self._scan: print self._scan
      for i, x in zip(self._crystal_ids, self._crystals): print i, x

      # Report on the refined parameters
      print self._param_report

    # Return the refinement history
    return self._refinery.history

  def selection_used_for_refinement(self):
    """Return a selection as a flex.bool in terms of the input reflection
    data of those reflections that were used in the final step of
    refinement."""

    from scitbx.array_family import flex
    matches = self._refman.get_matches()
    selection = flex.bool(self._refman.get_input_size(), False)

    try: # new reflection table format for matches
      isel = matches['iobs']
      selection.set_selected(isel, True)
    except TypeError: # old ObsPredMatch format for matches
      for m in matches:
        selection[m.iobs] = True

    return selection

  def predict_for_reflection_table(self, reflections):
    """perform prediction for all reflections in the supplied table"""

    # delegate to the target object, which has access to the predictor
    return self._target.predict_for_reflection_table(reflections)

  def predict_reflections(self):
    """Predict all reflection positions after refinement"""

    #FIXME
    raise NotImplementedError("predict_reflections is broken and to be deprecated")
    # code left here as a template for use of e.g. ScanVaryingReflectionListGenerator

    #FIXME only works for a single crystal
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.model.data import ReflectionList
    from math import pi
    from dials.algorithms.refinement.prediction.predictors import \
            ScanVaryingReflectionListGenerator

    if any([self._scan, self._goniometer]) is None:
        raise TypeError("Prediction can only be done when a scan and "
                        "goniometer are provided")

    s0 = self._beam.get_s0()
    dmin = self._detector.get_max_resolution(s0)

    # Test whether prediction is scan-varying or not (try-except duck-typing
    # fails because scan-varying prediction is done in multiple processes).
    from dials.algorithms.refinement.parameterisation import \
      VaryingCrystalPredictionParameterisation
    if isinstance(self._pred_param, VaryingCrystalPredictionParameterisation):

      sv_predictor = ScanVaryingReflectionListGenerator(self._crystal,
                            self._beam, self._goniometer, self._scan, dmin)
      refs = ReflectionList(sv_predictor())
      new_reflections = ray_intersection(self._detector, refs)

    else: # prediction seems to be scan-static
      from dials.algorithms.spot_prediction import IndexGenerator
      from cctbx.sgtbx import space_group_type
      from dials.algorithms.spot_prediction import RayPredictor
      from dials.algorithms.spot_prediction import reflection_frames

      # Create an index generator
      generate_hkl = IndexGenerator(self._crystal.get_unit_cell(),
                        space_group_type(self._crystal.get_space_group()),
                        dmin)

      # Create a spot predictor
      predict_rays = RayPredictor(s0,
                                  self._goniometer.get_rotation_axis(),
                                  self._scan.get_oscillation_range(deg=False))

      # predict
      miller_indices = generate_hkl.to_array()
      new_reflections = predict_rays(miller_indices, self._crystal.get_A())
      new_reflections = ray_intersection(self._detector, new_reflections)
      new_reflections = reflection_frames(self._scan, new_reflections)

    return new_reflections
