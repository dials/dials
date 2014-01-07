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
from dials.model.experiment.experiment_list import \
  ExperimentList, Experiment, ExperimentListFactory

class RefinerFactory(object):
  """Factory class to create refiners"""

  @staticmethod
  def from_parameters_data_experiments(params,
                                        reflections,
                                        experiments,
                                        verbosity=0):

    #TODO Checks on the input
    #E.g. does every experiment contain at least one overlapping model with at
    #least one other experiment? Are all the experiments either rotation series
    #or stills (the combination of both not yet supported)?

    # copy the experiments
    experiments = copy.deepcopy(experiments)

    # With this interface, assume that these come either from a scan, or
    # are None. The from_parameters_data_models interface provides more control.
    # We currently only support a single Experiment.
    scan = experiments[0].scan
    if scan:
      image_width_rad = scan.get_oscillation(deg=False)[1]
      sweep_range_rad = scan.get_oscillation_range(deg=False)
    else:
      image_width_rad = None
      sweep_range_rad = None

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
                params, experiments)

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
        experiments, sweep_range_rad, verbosity)

    return

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
                                  crystal_ids=None,
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
                params, experiments)

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
        experiments, sweep_range_rad, verbosity)

    if verbosity > 1:
      print ("Number of observations that pass initial inclusion criteria = %d"
             % refman.get_accepted_refs_size())
      print ("Working set size = %d observations"
             % refman.get_sample_size())
      print "Reflection manager built\n"

    if verbosity > 1: print "Building target function"

    # create target function
    target = RefinerFactory.config_target(params, crystals, crystal_ids, beam,
                    goniometer, detector, image_width_rad, refman,
                    pred_param)

    if verbosity > 1: print "Target function built\n"

    if verbosity > 1: print "Building refinement engine"

    # create refinery
    refinery = RefinerFactory.config_refinery(
                    params, target, pred_param, verbosity)

    if verbosity > 1: print "Refinement engine built\n"

    # build refiner interface and return
    return Refiner(reflections, beam, crystals, crystal_ids, detector,
                    pred_param, param_reporter, refman, target, refinery,
                    goniometer=goniometer,
                    scan=scan,
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

    # Shorten paths
    import dials.algorithms.refinement.parameterisation as par

    # FIXME: Multiple Experiments not yet supported!
    if len(experiments) > 1:
      raise RuntimeError("Multiple experiment parameterisation not"
                         "yet supported")
    beam = experiments[0].beam
    goniometer = experiments[0].goniometer
    crystal = experiments[0].crystal
    detector = experiments[0].detector
    scan = experiments[0].scan

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
    if crystal_options.scan_varying:
      if crystal_options.num_intervals == "fixed_width":
        sweep_range_deg = scan.get_oscillation_range(deg=True)
        deg_per_interval = crystal_options.interval_width_degrees
        n_intervals = int(
          abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval)
      else:
        n_intervals = crystal_options.absolute_num_intervals
      assert [goniometer, scan].count(None) == 0
      xl_ori_param = par.ScanVaryingCrystalOrientationParameterisation(
                                          crystal,
                                          scan.get_image_range(),
                                          n_intervals)
      xl_uc_param = par.ScanVaryingCrystalUnitCellParameterisation(
                                          crystal,
                                          scan.get_image_range(),
                                          n_intervals)
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
            track_parameter_correlation = options.track_parameter_correlation,
            max_iterations = options.max_iterations)

  @staticmethod
  def config_refman(params, reflections, experiments,
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

    #FIXME only single Experiment currently supported
    goniometer = experiments[0].goniometer
    beam = experiments[0].beam
    if goniometer:
      from dials.algorithms.refinement.target import ReflectionManager as refman

    else:
      from dials.algorithms.refinement.target_stills import \
          ReflectionManagerXY as refman

    # do outlier rejection?
    if options.do_outlier_rejection:
      iqr_multiplier=options.iqr_multiplier
    else:
      iqr_multiplier=None

    return refman(reflections=reflections,
                  beam=beam,
                  gonio=goniometer,
                  sweep_range_rad=sweep_range_rad,
                  nref_per_degree=nref_per_degree,
                  min_num_obs=options.minimum_number_of_reflections,
                  max_num_obs=options.maximum_number_of_reflections,
                  close_to_spindle_cutoff=options.close_to_spindle_cutoff,
                  iqr_multiplier=iqr_multiplier,
                  verbosity=verbosity)

  @staticmethod
  def config_target(params, crystals, crystal_ids, beam, goniometer, detector,
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
    if goniometer:
      from dials.algorithms.refinement.prediction import ReflectionPredictor
      ref_predictor = ReflectionPredictor(crystals, crystal_ids, beam, goniometer)

      import dials.algorithms.refinement.target as targ

      target = targ.LeastSquaresPositionalResidualWithRmsdCutoff(
                      ref_predictor, detector, refman, pred_param,
                      image_width_rad, options.bin_size_fraction,
                      absolute_cutoffs)
    else:
      from dials.algorithms.refinement.prediction import \
          StillsReflectionPredictor
      ref_predictor = StillsReflectionPredictor(crystals, crystal_ids, beam)

      import dials.algorithms.refinement.target_stills as targ

      target = targ.LeastSquaresXYResidualWithRmsdCutoff(
                      ref_predictor, detector, refman, pred_param,
                      options.bin_size_fraction,
                      absolute_cutoffs)

    return target

class Refiner(object):
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
    get_param_reporter
    parameter_correlation_plot
    selection_used_for_refinement
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

  def __init__(self, reflections, beam, crystals, crystal_ids, detector,
               pred_param, param_reporter, refman, target, refinery,
               goniometer=None,
               scan=None,
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

    # keep the data and models public for access after refinement
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

    self._refinery.run()

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
    selection = flex.bool(self._refman.get_input_size())
    for m in matches:
      selection[m.iobs] = True

    return selection

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
    from dials.algorithms.refinement.parameterisation import \
      VaryingCrystalPredictionParameterisation
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
