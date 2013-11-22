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

class RefinerFactory(object):
  """Factory class to create refiners"""

  @staticmethod
  def from_parameters_data_models(params,
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
                                  verbosity=0):
    """Given a set of parameters, reflections and experimental models, construct
    a refiner.

    Mandatory arguments:
      params - The input parameters as a phil scope_extract object
      reflections - Input ReflectionList data

    Argument alternatives:
      sweep - Experimental models as a dxtbx sweep
        or
      /beam - A dxtbx Beam object
      \detector - A dxtbx Detector object
      crystal - A dials Crystal object
        or
      crystals - A list of dials Crystal objects

    Optional arguments:
      goniometer - A dxtbx Goniometer model
      scan - A dxtbx Scan model
        or
      /image_width_rad - The 'width' of one image in radians
      \sweep_range_rad - Pair of rotation scan extremes in radians
      verbosity - An integer verbosity level

    Returns:
      The Refiner2 instance (to be renamed Refiner once the old
      interface is no longer in use).

    Notes
    -----

    The interface is intended to be flexible by allowing input to be passed in
    various forms. Some are incompatible, e.g. passing a sweep alongside any
    other dxtbx model is disallowed. Either crystal or crystals must be
    provided (but not both).

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
    if scan:
      image_width_rad = scan.get_oscillation(deg=False)[1]
      sweep_range_rad = scan.get_oscillation_range(deg=False)

    if goniometer and not image_width_rad:
      # if there is neither a scan nor the image width provided,
      # the target rmsd must be provided in absolute terms
      assert params.refinement.target.rmsd_cutoff == "absolute"
      # also there is no sweep range, so the reflection manager
      # cannot sample by number of reflections per degree
      assert params.refinement.reflections.reflections_per_degree is None

    # do we have the essential models?
    assert [beam, detector].count(None) == 0
    assert [crystal, crystals].count(None) == 1

    # copy the models
    from dxtbx.model import Beam, Detector
    from dxtbx.array_family import flex
    import copy
    # use copy constructors
    beam = Beam(beam)
    detector = Detector(flex.panel([panel for panel in detector]))
    if crystal: crystals = [copy.deepcopy(crystal)]
    if crystals: crystals = copy.deepcopy(crystals)

    # copy the reflections
    reflections = reflections.deep_copy()

    # check that the beam vectors are stored: if not, compute them
    from scitbx import matrix
    for ref in reflections:
      if ref.beam_vector != (0.0, 0.0, 0.0):
        continue
      panel = detector[ref.panel_number]
      x, y = panel.millimeter_to_pixel(ref.image_coord_mm)
      ref.beam_vector = matrix.col(panel.get_pixel_lab_coord(
          (x, y))).normalize() / beam.get_wavelength()

    # create parameterisations
    pred_param, param_reporter = \
            RefinerFactory.config_parameterisation(
                params, beam, detector, crystals, goniometer, scan)

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
    refman = RefinerFactory.config_refman(params, reflections,
        beam, goniometer, sweep_range_rad, verbosity)

    if verbosity > 1:
      print ("Number of observations that pass inclusion criteria = %d"
             % refman.get_accepted_refs_size())
      print ("Working set size = %d observations"
             % refman.get_sample_size())
      print "Reflection manager built\n"

    if verbosity > 1: print "Building target function"

    # create target function
    target = RefinerFactory.config_target(params, crystals, beam,
                    goniometer, detector, image_width_rad, refman,
                    pred_param)

    if verbosity > 1: print "Target function built\n"

    if verbosity > 1: print "Building refinement engine"

    # create refinery
    refinery = RefinerFactory.config_refinery(
                    params, target, pred_param, verbosity)

    if verbosity > 1: print "Refinement engine built\n"

    # build refiner interface and return
    return Refiner2(reflections, beam, crystals, detector,
                    pred_param, param_reporter, refman, target, refinery,
                    goniometer=goniometer,
                    scan=scan,
                    verbosity=verbosity)

  @staticmethod
  def config_parameterisation(
          params, beam, detector, crystals, goniometer, scan):
    """Given a set of parameters, create a parameterisation from a set of
    experimental models.

    Params:
        params The input parameters

    Returns:
        A tuple of the prediction equation parameterisation and the
        parameter reporter.
    """

    # Shorten parameter paths
    beam_options = params.refinement.parameterisation.beam
    crystal_options = params.refinement.parameterisation.crystal
    detector_options = params.refinement.parameterisation.detector
    prediction_options = params.refinement.parameterisation.prediction

    # Shorten paths
    import dials.algorithms.refinement.parameterisation as par

    # Beam (accepts goniometer=None)
    beam_param = par.BeamParameterisationOrientation(beam, goniometer)
    if beam_options.fix:
      if beam_options.fix == "all":
        beam_param.set_fixed([True, True])
      elif beam_options.fix == "in_spindle_plane":
        beam_param.set_fixed([True, False])
      else: # can only get here if refinement.phil is broken
        raise RuntimeError("beam_options.fix value not recognised")

    # Crystal
    # FIXME:
    if len(crystals) > 1:
      raise RuntimeError("Multiple crystal parameterisation not"
                         "yet supported")
    crystal = crystals[0] # Reminder for FIXME
    if crystal_options.scan_varying:
      assert [goniometer, scan].count(None) == 0
      xl_ori_param = par.ScanVaryingCrystalOrientationParameterisation(
                                          crystal,
                                          scan.get_image_range(),
                                          crystal_options.num_intervals)
      xl_uc_param = par.ScanVaryingCrystalUnitCellParameterisation(
                                          crystal,
                                          scan.get_image_range(),
                                          crystal_options.num_intervals)
    else:
      xl_ori_param = par.CrystalOrientationParameterisation(crystal)
      xl_uc_param = par.CrystalUnitCellParameterisation(crystal)

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

    # Detector
    if detector_options.panels == "automatic":
      if len(detector) > 1:
        det_param = par.DetectorParameterisationMultiPanel(detector, beam)
      else:
        det_param = par.DetectorParameterisationSinglePanel(detector)
    elif detector_options.panels == "single":
      det_param = DetectorParameterisationSinglePanel(detector)
    elif self._detector_par_options == "multiple":
      det_param = DetectorParameterisationMultiPanel(detector, beam)
    else: # can only get here if refinement.phil is broken
      raise RuntimeError("detector_options.panels value not recognised")

    if detector_options.fix:
      if detector_options.fix == "all":
        det_param.set_fixed([True] * det_param.num_total())
      elif detector_options.fix == "position":
        det_params = det_param.get_params(only_free = False)
        to_fix = [e.param_type.startswith('length') \
                  for e in det_params]
        det_param.set_fixed(to_fix)
      elif detector_options.fix == "orientation":
        det_params = det_param.get_params(only_free = False)
        to_fix = [e.param_type.startswith('angle') \
                  for e in det_params]
        det_param.set_fixed(to_fix)
      else: # can only get here if refinement.phil is broken
        raise RuntimeError("detector_options.fix value not recognised")

    # Prediction equation parameterisation
    crystal = crystals[0] # FIXME: multiple xls not yet supported
    if crystal_options.scan_varying:
      pred_param = par.VaryingCrystalPredictionParameterisation(
          detector, beam, crystal, goniometer,
          [det_param], [beam_param], [xl_ori_param], [xl_uc_param])
    elif goniometer is None:
      pred_param = par.XYPredictionParameterisation(
          detector, beam, crystal, goniometer,
          [det_param], [beam_param], [xl_ori_param], [xl_uc_param])
    else:
      pred_param = par.XYPhiPredictionParameterisation(
          detector, beam, crystal, goniometer,
          [det_param], [beam_param], [xl_ori_param], [xl_uc_param])

    # Parameter reporting
    param_reporter = par.ParameterReporter([det_param], [beam_param],
        [xl_ori_param], [xl_uc_param])

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

    if options.engine  == "SimpleLBFGS":
      from engine import SimpleLBFGS as refinery
    elif options.engine == "LBFGScurvs":
      from engine import LBFGScurvs as refinery
    elif options.engine == "GaussNewtonIterations":
      from engine import GaussNewtonIterations as refinery
    elif options.engine == "LevMarIterations":
      from engine import LevenbergMarquardtIterations as refinery
    else:
      raise RuntimeError("Refinement engine " + options.engine +
                         " not recognised")

    return refinery(target = target,
                    prediction_parameterisation = pred_param,
                    log = options.log,
                    verbosity = verbosity,
                    track_step = options.track_step,
                    track_gradient = options.track_gradient,
                    max_iterations = options.max_iterations)

  @staticmethod
  def config_refman(params, reflections, beam, goniometer,
                       sweep_range_rad, verbosity):
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
      if verbosity > 1:
        print "Random seed set to %d\n" % options.random_seed

    if goniometer:
      from dials.algorithms.refinement.target import ReflectionManager as refman

    else:
      from dials.algorithms.refinement.target_stills import \
          ReflectionManagerXY as refman

    return refman(reflections=reflections,
                  beam=beam,
                  gonio=goniometer,
                  sweep_range_rad=sweep_range_rad,
                  nref_per_degree=nref_per_degree,
                  min_num_obs=options.minimum_number_of_reflections,
                  max_num_obs=options.maximum_number_of_reflections,
                  inclusion_cutoff=options.inclusion_cutoff,
                  verbosity=verbosity)

  @staticmethod
  def config_target(params, crystals, beam, goniometer, detector,
      image_width_rad, refman, pred_param):
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

    # Determine whether the target is in X, Y, Phi space or just X, Y.
    crystal = crystals[0] # FIXME: multiple crystals not yet supported
    if goniometer:
      from dials.algorithms.refinement.prediction import ReflectionPredictor
      ref_predictor = ReflectionPredictor(crystal, beam, goniometer)

      import dials.algorithms.refinement.target as targ

      target = targ.LeastSquaresPositionalResidualWithRmsdCutoff(
                      ref_predictor, detector, refman, pred_param,
                      image_width_rad, options.bin_size_fraction,
                      absolute_cutoffs)
    else:
      from dials.algorithms.refinement.prediction import \
          StillsReflectionPredictor
      ref_predictor = StillsReflectionPredictor(crystal, beam)

      import dials.algorithms.refinement.target_stills as targ

      target = targ.LeastSquaresXYResidualWithRmsdCutoff(
                      ref_predictor, detector, refman, pred_param,
                      options.bin_size_fraction,
                      absolute_cutoffs)

    return target

  # FIXME Deprecated method used by old interface only
  @staticmethod
  def from_parameters(params, verbosity):
    """Given a set of parameters, construct the refiner

    Params:
        params The input parameters
        verbosity The verbosity level

    Returns:
        The refiner instance

    """

    parameterisation_strategy = \
                RefinerFactory.configure_parameterisation(params)
    refinery_strategy = RefinerFactory.configure_refinery(params)
    reflections_strategy = RefinerFactory.configure_refman(params)
    target_strategy = RefinerFactory.configure_target(params)

    return Refiner(parameterisation_strategy, refinery_strategy,
             reflections_strategy, target_strategy, verbosity)

  # FIXME Deprecated method used by old interface only
  @staticmethod
  def configure_parameterisation(params):
    """Given a set of parameters, configure a factory to build a
    parameterisation from a set of experimental models

    Params:
        params The input parameters

    Returns:
        The parameterisation factory instance
    """

    # Shorten parameter paths
    beam_options = params.refinement.parameterisation.beam
    crystal_options = params.refinement.parameterisation.crystal
    detector_options = params.refinement.parameterisation.detector
    prediction_options = params.refinement.parameterisation.prediction

    return ParameterisationFactory(beam_options, crystal_options,
                            detector_options, prediction_options)

  # FIXME Deprecated method used by old interface only
  @staticmethod
  def configure_refinery(params):
    """Given a set of parameters, configure a factory to build a
    refinery

    Params:
        params The input parameters

    Returns:
        The refinery factory instance
    """

    # Shorten parameter path
    options = params.refinement.refinery
    return RefineryFactory(options)

  # FIXME Deprecated method used by old interface only
  @staticmethod
  def configure_refman(params):
    """Given a set of parameters, configure a factory to build a
    reflection manager

    Params:
        params The input parameters

    Returns:
        The reflection manager factory instance
    """

    # Shorten parameter path
    options = params.refinement.reflections
    return RefmanFactory(options)

  # FIXME Deprecated method used by old interface only
  @staticmethod
  def configure_target(params):
    """Given a set of parameters, configure a factory to build a
    target function

    Params:
        params The input parameters

    Returns:
        The target factory instance
    """

    # Shorten parameter path
    options = params.refinement.target
    return TargetFactory(options)

# FIXME Deprecated method used by old interface only
class ParameterisationFactory(object):
  """ Factory class to create beam, crystal and detector parameterisations
  plus a parameterisation of the prediction equation."""

  def __init__(self, beam_options, crystal_options, detector_options,
               prediction_options):

    # Shorten paths
    import dials.algorithms.refinement.parameterisation as par

    # Beam
    self._beam_par = par.BeamParameterisationOrientation
    self._beam_fix = beam_options.fix

    # Crystal
    self._crystal_fix = crystal_options.fix
    self._crystal_scan_varying = crystal_options.scan_varying
    self._crystal_num_intervals = crystal_options.num_intervals

    if self._crystal_scan_varying:
      cop = par.ScanVaryingCrystalOrientationParameterisation
      cucp = par.ScanVaryingCrystalUnitCellParameterisation
    else:
      cop = par.CrystalOrientationParameterisation
      cucp = par.CrystalUnitCellParameterisation
    self._crystal_ori_par = cop
    self._crystal_uc_par = cucp

    # Detector
    if detector_options.panels not in ["automatic", "single", "multiple"]:
      raise RuntimeError("detector parameterisation type not recognised")

    self._detector_par_options = detector_options.panels
    self._detector_fix = detector_options.fix

    # Prediction equation parameterisation
    self._prediction_par_options = prediction_options.space
    if self._crystal_scan_varying:
      pep = par.VaryingCrystalPredictionParameterisation
    elif self._prediction_par_options == "XYPhi":
      pep = par.XYPhiPredictionParameterisation
    elif self._prediction_par_options == "XY":
      pep = par.XYPredictionParameterisation
    else:
      raise RuntimeError("Prediction equation type " +
          self._prediction_par_options + " not recognised")

    self.prediction_par = pep

    # Parameter reporting
    self.param_reporter = par.ParameterReporter

  def __call__(self, beam, crystal, goniometer, detector, scan):

    beam_param = self._beam_par(beam, goniometer)
    if self._beam_fix:
      if self._beam_fix == "all":
        beam_param.set_fixed([True, True])
      elif self._beam_fix == "in_spindle_plane":
        beam_param.set_fixed([True, False])

    if self._crystal_scan_varying:
      xl_ori_param = self._crystal_ori_par(crystal,
                                           scan.get_image_range(),
                                           self._crystal_num_intervals)
      xl_uc_param = self._crystal_uc_par(crystal,
                                         scan.get_image_range(),
                                         self._crystal_num_intervals)
    else:
      xl_ori_param = self._crystal_ori_par(crystal)
      xl_uc_param = self._crystal_uc_par(crystal)

    if self._crystal_fix:
      if self._crystal_fix == "all":
        xl_ori_param.set_fixed([True] * xl_ori_param.num_total())
        xl_uc_param.set_fixed([True] * xl_uc_param.num_total())
      elif self._crystal_fix == "cell":
        xl_uc_param.set_fixed([True] * xl_uc_param.num_total())
      elif self._crystal_fix == "orientation":
        xl_ori_param.set_fixed([True] * xl_ori_param.num_total())

    from dials.algorithms.refinement.parameterisation.detector_parameters \
        import DetectorParameterisationSinglePanel, \
            DetectorParameterisationMultiPanel

    if self._detector_par_options == "automatic":
      if len(detector) > 1:
        det_param = DetectorParameterisationMultiPanel(detector, beam)
      else:
        det_param = DetectorParameterisationSinglePanel(detector)
    if self._detector_par_options == "single":
      det_param = DetectorParameterisationSinglePanel(detector)
    if self._detector_par_options == "multiple":
      det_param = DetectorParameterisationMultiPanel(detector, beam)

    if self._detector_fix:
      if self._detector_fix == "all":
        det_param.set_fixed([True] * det_param.num_total())
      elif self._detector_fix == "position":
        det_params = det_param.get_params(only_free = False)
        to_fix = [e.param_type.startswith('length') \
                  for e in det_params]
        det_param.set_fixed(to_fix)
      elif self._detector_fix == "orientation":
        det_params = det_param.get_params(only_free = False)
        to_fix = [e.param_type.startswith('angle') \
                  for e in det_params]
        det_param.set_fixed(to_fix)

    pred_param = self.prediction_par(detector, beam, crystal, goniometer,
            [det_param], [beam_param], [xl_ori_param], [xl_uc_param])

    param_reporter = self.param_reporter([det_param], [beam_param],
        [xl_ori_param], [xl_uc_param])

    return (beam_param, xl_ori_param, xl_uc_param, det_param, pred_param,
            param_reporter)

# FIXME Deprecated method used by old interface only
class RefineryFactory(object):
  """Factory class to create a Refinery object (the refinement engine)"""

  def __init__(self, options):

    import dials.algorithms.refinement.engine as engine

    if options.engine  == "SimpleLBFGS":
      from engine import SimpleLBFGS as ref
    elif options.engine == "LBFGScurvs":
      from engine import LBFGScurvs as ref
    elif options.engine == "GaussNewtonIterations":
      from engine import GaussNewtonIterations as ref
    elif options.engine == "LevMarIterations":
      from engine import LevenbergMarquardtIterations as ref
    else:
      raise RuntimeError("Refinement engine " + options.engine +
                         " not recognised")

    self._refinery = ref
    self._track_step = options.track_step
    self._track_gradient = options.track_gradient
    self._logfile = options.log
    self._max_iterations = options.max_iterations

  def __call__(self, target, prediction_parameterisation, verbosity):

    return self._refinery(
        target = target,
        prediction_parameterisation = prediction_parameterisation,
        log = self._logfile,
        verbosity = verbosity,
        track_step = self._track_step,
        track_gradient = self._track_gradient,
        max_iterations = self._max_iterations)

# FIXME Deprecated method used by old interface only
class RefmanFactory(object):
  """Factory class to create a ReflectionManager"""

  def __init__(self, options):

    if options.implementation == "rotation":
      from dials.algorithms.refinement.target import \
        ReflectionManager as refman
    elif options.implementation == "stills":
      from dials.algorithms.refinement.target_stills import \
        ReflectionManagerXY as refman
    else:
      raise RuntimeError("ReflectionManager type " +
                         options.implementation + " not recognised")
    self._refman = refman

    self._random_seed = options.random_seed

    self._ref_per_degree = options.reflections_per_degree
    if options.use_all_reflections:
      self._ref_per_degree = None

    self._max_num_obs = options.maximum_number_of_reflections

    self._min_num_obs = options.minimum_number_of_reflections

    self._inclusion_cutoff = options.inclusion_cutoff

  def __call__(self, reflections, beam, goniometer, scan, verbosity):

    # While a random subset of reflections is used, continue to
    # set random.seed to get consistent behaviour
    if self._random_seed:
      import random
      random.seed(self._random_seed)
      if verbosity > 1:
        print "Random seed set to %d\n" % self._random_seed

    sweep_range = scan.get_oscillation_range(deg=False)
    return self._refman(reflections=reflections,
                        beam=beam,
                        gonio=goniometer,
                        sweep_range_rad=sweep_range,
                        nref_per_degree=self._ref_per_degree,
                        min_num_obs=self._min_num_obs,
                        max_num_obs=self._max_num_obs,
                        inclusion_cutoff=self._inclusion_cutoff,
                        verbosity=verbosity)

# FIXME Deprecated method used by old interface only
class TargetFactory(object):
  """Factory class to create a target function object"""

  def __init__(self, options):

    # this is a temporary patch for the old interface
    self._XY = False

    if options.implementation == "basic":
      from dials.algorithms.refinement.target import \
          LeastSquaresPositionalResidualWithRmsdCutoff as targ
    elif options.implementation == "XY":
      from dials.algorithms.refinement.target_stills import \
          LeastSquaresXYResidualWithRmsdCutoff as targ
      self._XY = True
    else:
      raise RuntimeError("Target type " + options.implementation +
                          " not recognised")

    self._frac_binsize_cutoff = options.bin_size_fraction
    if options.rmsd_cutoff == "fraction_of_bin_size":
      self._absolute_cutoffs = None
    elif options.rmsd_cutoff == "absolute":
      self._absolute_cutoffs = options.absolute_cutoffs
    else:
      raise RuntimeError("Target function rmsd_cutoff option" +
          options.rmsd_cutoff + " not recognised")

    # Reflection prediction
    from dials.algorithms.refinement.prediction import \
        ReflectionPredictor as rp

    self._target = targ
    self._ref_predictor = rp

  def __call__(self, crystal, beam, goniometer, detector, scan,
      refman, pred_param):

    image_width = scan.get_oscillation(deg=False)[1]

    rp = self._ref_predictor(crystal, beam, goniometer)
    if self._XY:
      return self._target(rp, detector, refman, pred_param,
                      self._frac_binsize_cutoff,
                      self._absolute_cutoffs)
    else:
      return self._target(rp, detector, refman, pred_param,
                      image_width, self._frac_binsize_cutoff,
                      self._absolute_cutoffs)

class Refiner2(object):
  """Public interface for performing DIALS refinement.

  Public methods:
    run
    rmsds
    get_beam
    get_crystal
    get_crystals
    get_detector
    get_goniometer
    get_scan
    get_reflections
    get_matches
    selection_used_for_refinement
    write_residuals_table
    predict_reflections

  Notes:
    * The return value of run is a recorded history of the refinement
    * The model accessors provide copies of those models that might be modified
      by refinement (beam, crystal and detector)
    * get_matches exposes the function of the same name from the privately
      stored reflection manager
    * The return value of selection_used_for_refinement is a flex.bool
    * predict_reflections does so for all observable reflections, not just the
      set passed into the refiner

    """

  def __init__(self, reflections, beam, crystals, detector,
               pred_param, param_reporter, refman, target, refinery,
               goniometer=None,
               scan=None,
               verbosity=0):
    """
    Mandatory arguments:
      reflections - Input ReflectionList data
      beam - A dxtbx Beam object
      crystals - A list of DIALS Crystal objects
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

    # keep the data and models public for access after refinement
    self._reflections = reflections
    self._beam = beam
    self._crystals = crystals
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
    from dxtbx.array_family import flex
    return Detector(flex.panel([deepcopy(panel) \
                                for panel in self._detector]))

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
    """Delegated to the reflection manager, but suppress output"""

    # FIXME Consider: Does this information really need to be exposed by the
    # public API (indexing code seems to use it, but is it necessary?)
    return self._refman.get_matches(silent = True)

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
      for x in self._crystals: print x

    self._refinery.run()

    if self._verbosity > 1:
      print
      print "Experimental models after refinement"
      print "------------------------------------"
      print self._beam
      print self._detector
      if self._goniometer: print self._goniometer
      if self._scan: print self._scan
      for x in self._crystals: print x

      # Write scan-varying parameters to file, if there were any
      if self._scan and self._param_report.varying_params_vs_image_number(
              self._scan.get_image_range()):
        print "Writing scan-varying parameter table to file"

      # Report on the refined parameters
      print self._param_report

      print "Writing residuals to file"
      self.write_residuals_table()

    # Return the refinement history
    return self._refinery.history

  def selection_used_for_refinement(self):
    """Return a selection as a flex.bool in terms of the input reflection
    data of those reflections that were used in the final step of
    refinement."""

    from scitbx.array_family import flex
    matches = self._refman.get_matches()
    selection = flex.bool(self._refman.get_input_size())
    for m in matches:
      selection[m.iobs] = True

    return selection

  def write_residuals_table(self):

    matches = self._refman.get_matches()

    f = open("residuals.dat","w")
    header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
        "Y_calc\tPhi_calc\n")
    f.write(header)

    for m in matches:
      msg = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%9.6f\t"
            "%5.3f\n")
      msg = msg % (m.H[0], m.H[1], m.H[2], m.frame_o, m.Xo, m.Yo,
                   m.Phio, m.Xc, m.Yc, m.Phic)
      f.write(msg)
    f.close()

  def predict_reflections(self):
    """Predict all reflection positions after refinement"""

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
    if isinstance(self._pred_param, VaryingCrystalPredictionParameterisation):

      sv_predictor = ScanVaryingReflectionListGenerator(self._pred_param,
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

class Refiner(object):
  """The refiner class."""

  def __init__(self, parameterisation_strategy, refinery_strategy,
               reflections_strategy, target_strategy, verbosity):
    """ Initialise the refiner class.

    Params:
        parameterisation_strategy The parameterisation strategy
        refinery_strategy The engine strategy
        reflections_strategy The reflection manager strategy
        target_strategy The target function strategy
        verbosity The verbosity level
    """

    self.create_param = parameterisation_strategy
    self.create_refinery = refinery_strategy
    self.create_target = target_strategy
    self.create_refman = reflections_strategy
    self._verbosity = verbosity

  def prepare(self, sweep, crystal, reflections):
    """ Prepare refiner with experimental models and data.

    Params:
        sweep The sweep to process
        crystal The crystal to process
        reflections The reflection list

    Returns:
        The refined (sweep, crystal)

    """
    from scitbx import matrix

    # Get the models from the sweep
    self.sweep = sweep
    self.beam = sweep.get_beam()
    self.detector = sweep.get_detector()
    self.gonio = sweep.get_goniometer()
    self.scan = sweep.get_scan()
    self.crystal = crystal

    if self._verbosity > 1:
      print ""
      print "Experimental Models"
      print "-------------------"
      print self.beam
      print self.detector
      print self.gonio
      print self.scan
      print self.crystal


    # Copy the reflections
    self.reflections = reflections
    self._saved_reflections = self.reflections.deep_copy()

    # check that the beam vectors are stored: if not, compute them
    for ref in self.reflections:
      if ref.beam_vector != (0.0, 0.0, 0.0):
        continue
      panel = self.detector[ref.panel_number]
      x, y = panel.millimeter_to_pixel(ref.image_coord_mm)
      ref.beam_vector = matrix.col(panel.get_pixel_lab_coord(
          (x, y))).normalize() / self.beam.get_wavelength()

    ###########################
    # Parameterise the models #
    ###########################

    (self.beam_param, self.xl_ori_param, self.xl_uc_param, self.det_param,
     self.pred_param, self.param_report) = \
        self.create_param(self.beam, self.crystal, self.gonio,
                          self.detector, self.scan)

    if self._verbosity > 1:
      print "Prediction equation parameterisation built\n"
      print "Parameter order : name mapping"
      for i, e in enumerate(self.pred_param.get_param_names()):
        print "Parameter %03d : " % i + e
      print

      print "Prior to refinement the experimental model is:"
      print_model_geometry(self.beam, self.detector, self.crystal)
      print

    #####################################
    # Select reflections for refinement #
    #####################################

    if self._verbosity > 1:
      print "Building reflection manager"
      print ("Input reflection list size = %d observations"
             % len(self.reflections))

    self.refman = self.create_refman(self.reflections, self.beam,
                            self.gonio, self.scan, self._verbosity)

    if self._verbosity > 1:
      print ("Number of observations that pass inclusion criteria = %d"
             % self.refman.get_accepted_refs_size())
      print ("Working set size = %d observations"
             % self.refman.get_sample_size())
      print "Reflection manager built\n"

    ##############################
    # Set up the target function #
    ##############################

    if self._verbosity > 1: print "Building target function"

    self.target = self.create_target(self.crystal, self.beam,
        self.gonio, self.detector, self.scan, self.refman,
        self.pred_param)

    if self._verbosity > 1: print "Target function built\n"

    ################################
    # Set up the refinement engine #
    ################################

    if self._verbosity > 1: print "Building refinement engine"

    self.refinery = self.create_refinery(self.target, self.pred_param,
                                         self._verbosity)

    if self._verbosity > 1: print "Refinement engine built\n"

    return

  def rmsds(self):
    """Return rmsds of the current model"""

    self.refinery.prepare_for_step()

    return self.target.rmsds()

  def __call__(self, sweep=None, crystal=None, reflections=None):
    """Run refinement"""

    if sweep and crystal and reflections:
      self.prepare(sweep, crystal, reflections)
    else: assert [sweep, crystal, reflections].count(None) == 3

    ###################################
    # Do refinement and return models #
    ###################################

    self.refinery.run()

    if self._verbosity > 1:
      print
      print "Refinement has completed with the following geometry:"
      print_model_geometry(self.beam, self.detector, self.crystal)

      if self.param_report.varying_params_vs_image_number(
          self.scan.get_image_range()):
        print "Writing scan-varying parameter table to file"

      print "Reporting on the refined parameters:"
      print self.param_report

      print "Writing residuals to file"
      self.write_residuals_table()

    # Do a test of new reflection pos
    #self._update_reflections_test()

    # Return the refinery, containing useful information such as the
    # refinement history. The refined models are set by side-effect
    return self.refinery

  def selection_used_for_refinement(self):
    """Return a selection as a flex.bool in terms of the input reflection
    data of those reflections that were used in the final step of
    refinement."""

    from scitbx.array_family import flex
    matches = self.refman.get_matches()
    selection = flex.bool(self.refman.get_input_size())
    for m in matches:
      selection[m.iobs] = True

    return selection

  def write_residuals_table(self):

    matches = self.refman.get_matches()

    f = open("residuals.dat","w")
    header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
        "Y_calc\tPhi_calc\n")
    f.write(header)

    for m in matches:
      msg = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%9.6f\t"
            "%5.3f\n")
      msg = msg % (m.H[0], m.H[1], m.H[2], m.frame_o, m.Xo, m.Yo,
                   m.Phio, m.Xc, m.Yc, m.Phic)
      f.write(msg)
    f.close()

  def _update_reflections_test(self):
    from cctbx.array_family import flex
    from collections import defaultdict

    # Get miller indices from saved reflectons
    miller_indices = [r.miller_index for r in self._saved_reflections]

    self.miller_indices = flex.miller_index(miller_indices)

    print "Predicting new reflections"
    self.predict_reflections()

    # Put coords from same hkl in dict for saved reflections
    coord1 = defaultdict(list)
    for r1 in self._saved_reflections:
      coord1[r1.miller_index].append(r1.image_coord_px)

    # Put coords from same hkl in dict for new reflections
    coord2 = defaultdict(list)
    for r2 in self._new_reflections:
      coord2[r2.miller_index].append(r2.image_coord_px)

    # Print out coords for each hkl
    for h in coord1.keys():
      c1 = coord1[h]
      c2 = coord2[h]
      #print c1, c2

  def predict_reflections(self):
    """Predict all reflection positions after refinement and make the
    bounding boxes."""
    from dials.algorithms.integration import ReflectionPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.model.data import ReflectionList
    from math import pi
    from dials.algorithms.refinement.prediction.predictors import \
            ScanVaryingReflectionListGenerator

    s0 = self.beam.get_s0()
    dmin = self.detector.get_max_resolution(s0)
    sv_predictor = ScanVaryingReflectionListGenerator(self.pred_param,
                        self.beam, self.gonio, self.scan, dmin)

    # Duck typing to determine whether prediction is scan-varying or not
    try:
      refs = ReflectionList(sv_predictor())
      self._new_reflections = ray_intersection(self.detector, refs)

    except AttributeError: # prediction seems to be scan-static
      predict = ReflectionPredictor()
      self._new_reflections = predict(self.sweep, self.crystal)

    self.sigma_divergence = self.beam.get_sigma_divergence()
    self.sigma_mosaicity = self.crystal.get_mosaicity()

    # Set the divergence and mosaicity
    n_sigma = 5.0
    delta_divergence = n_sigma * self.sigma_divergence * pi / 180.0
    delta_mosaicity = n_sigma * self.sigma_mosaicity * pi / 180.0

    # FIXME: DIALS_ASSERT(delta_divergence > 0.0) failure.
    #
    # Create the bounding box calculator
    #calculate_bbox = BBoxCalculator(self.beam, self.detector, self.gonio,
    #    self.scan, delta_divergence, delta_mosaicity)

    # Calculate the frame numbers of all the reflections
    #calculate_bbox(self._new_reflections)

    return self._new_reflections
