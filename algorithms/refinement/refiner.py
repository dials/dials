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
from logging import info, debug, warning

from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
from dials.array_family import flex
from libtbx.phil import parse
from libtbx.utils import Sorry
import libtbx

# Include external phil scopes as strings to avoid problems
# with the include scope directive
from dials.algorithms.refinement.outlier_detection.outlier_base \
  import phil_str as outlier_phil_str
from dials.algorithms.refinement.restraints.restraints_parameterisation \
  import uc_phil_str as uc_restraints_phil_str
format_data = {'outlier_phil':outlier_phil_str, 'uc_restraints_phil':uc_restraints_phil_str}
phil_scope = parse('''

refinement
  .help = "Parameters to configure the refinement"
{

  mp
    .expert_level = 1
  {
    nproc = 1
      .type = int(value_min=1)
      .help = "The number of processes to use. Only applicable to certain"
              "choices of refinement engine!"
  }

  verbosity = 1
    .help = "verbosity level"
    .type = int(value_min=0)

  parameterisation
    .help = "Parameters to control the parameterisation of experimental models"
  {

    auto_reduction
      .help = "determine behaviour when there are too few reflections to"
              "reasonably produce a full parameterisation of the experiment list"
    {
      min_nref_per_parameter = 5
        .help = "the smallest number of reflections per parameter for a"
                "model parameterisation below which the parameterisation will"
                "not be made in full, but the action described below will be"
                "triggered."
        .type = int(value_min=1)

      action = *fail fix remove
        .help = "action to take if there are too few reflections across the"
                "experiments related to a particular model parameterisation."
                "If fail, an exception will be raised and refinement will not"
                "proceed. If fix, refinement will continue but with the"
                "parameters relating to that model remaining fixed at their"
                "initial values. If remove, parameters relating to that model"
                "will be fixed, and in addition all reflections related to"
                "that parameterisation will be removed. This will therefore"
                "remove this reflections from other parameterisations of the"
                "global model too. For example, if a crystal model could not"
                "be parameterised it will be excised completely and not"
                "contribute to the joint refinement of the detector and beam."
                "In the fix mode, reflections emanating from that crystal will"
                "still form residuals and will contribute to detector and beam"
                "refinement."
        .type = choice
    }

    beam
      .help = "beam parameters"
    {
      fix = all *in_spindle_plane out_spindle_plane *wavelength
        .help = "Whether to fix beam parameters. By default,"
                "in_spindle_plane is selected, and one of the two"
                "parameters is fixed. If a goniometer is present"
                "this leads to the beam orientation being restricted"
                "to a direction in the initial spindle-beam plane."
                "Wavelength is also fixed by default, to allow refinement of"
                "the unit cell volume."
        .type = choice(multi=True)

      fix_list = None
        .type = ints(value_min=0)
        .help = "Fix specified parameters by a list of indices"
        .expert_level = 1
    }

    crystal
      .help = "crystal parameters"
    {
      fix = all cell orientation
        .help = "Fix crystal parameters"
        .type = choice

      unit_cell
      {
        fix_list = None
          .type = ints(value_min=0)
          .help = "Fix specified parameters by a list of indices"
          .expert_level = 1

        %(uc_restraints_phil)s
      }

      orientation
      {
        fix_list = None
          .type = ints(value_min=0)
          .help = "Fix specified parameters by a list of indices"
          .expert_level = 1
      }

      scan_varying = False
        .help = "Parameterise the crystal to vary during the scan"
        .type = bool

      num_intervals = *fixed_width absolute
        .help = "Choose the way to determine the number of intervals for scan-"
                "varying refinement"
        .type = choice
        .expert_level = 1

      interval_width_degrees = 36.0
        .help = "Width of scan between checkpoints in degrees"
        .type = float(value_min=0.)
        .expert_level = 1

      absolute_num_intervals = 5
        .help = "Number of intervals between checkpoints if scan_varying"
                "refinement is requested"
        .type = int(value_min=1)
        .expert_level = 1

      UB_model_per = reflection image *block
        .help = "Compose a new crystal model either every reflection (slow),"
                "every image (faster, less accurate) or within blocks of a"
                "width specified in the reflections parameters. When this block"
                "width is larger than the image width the result is faster"
                "again, with another trade-off in accuracy"
        .type = choice
        .expert_level = 1
     }

    detector
      .help = "detector parameters"
    {
      panels = *automatic single multiple hierarchical
        .help = "Select appropriate detector parameterisation. Both the"
                "single and multiple panel detector options treat the whole"
                "detector as a rigid body. The hierarchical parameterisation"
                "treats groups of panels as separate rigid bodies."
        .type = choice
        .expert_level = 1

      hierarchy_level = 0
        .help = "Level of the detector hierarchy (starting from the root at 0)"
                "at which to determine panel groups to parameterise independently"
        .type = int(value_min=0)
        .expert_level = 1

      fix = all position orientation
        .help = "Fix detector parameters. The translational parameters"
                "(position) may be set separately to the orientation."
        .type = choice

      fix_list = None
        .type = ints(value_min=0)
        .help = "Fix specified parameters by a list of indices"
        .expert_level = 1
    }

    sparse = Auto
      .help = "Calculate gradients using sparse data structures."
      .type = bool
      .expert_level = 1

    treat_single_image_as_still = False
      .help = "Set this to True to treat a single image scan with a non zero"
              "oscillation width as a still"
      .type = bool
      .expert_level = 1
  }

  refinery
    .help = "Parameters to configure the refinery"
    .expert_level = 1
  {
    engine = SimpleLBFGS LBFGScurvs GaussNewton *LevMar SparseLevMar
      .help = "The minimisation engine to use"
      .type = choice

    track_step = False
      .help = "Record parameter shifts history in the refinement journal, if"
              "the engine supports it."
      .type = bool

    track_gradient = False
      .help = "Record parameter gradients history in the refinement journal, if"
              "the engine supports it."
      .type = bool

    track_parameter_correlation = False
      .help = "Record correlation matrix between columns of the Jacobian for"
              "each step of refinement."
      .type = bool

    track_out_of_sample_rmsd = False
      .type = bool
      .help = "Record RMSDs calculated using the refined experiments with"
              "reflections not used in refinement at each step. Only valid if a"
              "subset of input reflections was taken for refinement"

    log = None
      .help = "Filename for an optional log that a minimisation engine may use"
              "to write additional information"
      .type = path

    max_iterations = None
      .help = "Maximum number of iterations in refinement before termination."
              "None implies the engine supplies its own default."
      .type = int(value_min=1)
  }

  target
    .help = "Parameters to configure the target function"
    .expert_level = 1
  {

    rmsd_cutoff = *fraction_of_bin_size absolute
      .help = "Method to choose rmsd cutoffs. This is currently either as a"
              "fraction of the discrete units of the spot positional data, i.e."
              "(pixel width, pixel height, image thickness in phi), or a tuple"
              "of absolute values to use as the cutoffs"
      .type = choice

    bin_size_fraction = 0.2
      .help = "Cut off in the natural discrete units of positional data, viz.,"
              "(pixel width, pixel height, image thickness in phi) to use to"
              "determine when the RMSD target is achieved. Only used if"
              "rmsd_cutoff = fraction_of_bin_size."
      .type = float(value_min=0.)

    absolute_cutoffs = None
      .help = "Absolute Values for the RMSD target achieved cutoffs in X, Y and"
              "Phi. The units are (mm, mm, rad)."
      .type = floats(size=3, value_min=0.)

    gradient_calculation_blocksize = None
      .help = "Maximum number of reflections to use for gradient calculation."
              "If there are more reflections than this in the manager then"
              "the minimiser must do the full calculation in blocks."
      .type = int(value_min=1)

  }

  reflections
    .help = "Parameters used by the reflection manager"
  {

    reflections_per_degree = 100
      .help = "The number of centroids per degree of the sweep to use in"
              "refinement."
      .type = float(value_min=0.)

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input"
              "reflection list is done"
      .type = int

    maximum_sample_size = None
      .help = "The maximum number of reflections to use in refinement."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)

    use_all_reflections = False
      .help = "Override reflections_per_degree and use all available centroids"
              "in refinement."
      .type = bool

    random_seed = 42
      .help = "Random seed to use when sampling to create a working set of"
              "reflections. May be int or None."
      .type = int
      .expert_level = 1

    close_to_spindle_cutoff = 0.02
      .help = "The inclusion criterion currently uses the volume of the"
              "parallelepiped formed by the spindle axis, the incident"
              "beam and the scattered beam. If this is lower than some"
              "value then the reflection is excluded from refinement."
              "In detector space, these are the reflections located close"
              "to the rotation axis."
      .type = float(value_min = 0)
      .expert_level = 1

    block_width = 1.0
      .help = "Width of a reflection 'block' (in degrees) determining how fine-"
              "grained the model used for scan-varying prediction during"
              "refinement is. Currently only has any effect if the crystal"
              "parameterisation is set to use UB_model_per=block"
      .type = float(value_min = 0.)
      .expert_level = 1

    weighting_strategy
      .help = "Parameters to configure weighting strategy overrides"
      .expert_level = 1
    {
      override = statistical stills constant
        .help = "selection of a strategy to override default weighting behaviour"
        .type = choice

      delpsi_constant = 10000
        .help = "used by the stills strategy to choose absolute weight value"
                "for the angular distance from Ewald sphere term of the target"
                "function, whilst the X and Y parts use statistical weights"
        .type = float(value_min = 0)

      constants = 1.0 1.0 1.0
        .help = "constant weights for three parts of the target function,"
                "whether the case is for stills or scans. The default gives"
                "unit weighting."
        .type = floats(size = 3, value_min = 0)
    }

    # Does not work:
    #
    # include scope dials.algorithms.refinement.outlier_detection.phil_scope
    #
    # instead just paste the string in directly here.
    %(outlier_phil)s
  }
}
'''%format_data, process_includes=True)

class RefinerFactory(object):
  """Factory class to create refiners"""

  @classmethod
  def _filter_reflections(cls, reflections):
    '''Return a copy of the input reflections filtered to keep only
    those columns that are required by refinement'''

    cols = ['id', 'miller_index', 'panel', 's1', 'xyzobs.mm.value',
            "xyzobs.px.value", "xyzcal.px",
            'xyzobs.mm.variance', 'flags']
    # NB xyzobs.px.value & xyzcal.px required by SauterPoon outlier rejector
    rt = flex.reflection_table()

    # copy columns to the new table. Could use the select method
    # for this except that 's1' is optional in the input so would want
    # to copy that in like this if present anyway
    for k in cols:
      if k in reflections.keys():
        rt[k] = reflections[k]
    return rt

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
    if verbosity is None:
      verbosity = params.refinement.verbosity

    # copy the experiments
    import copy
    experiments = copy.deepcopy(experiments)

    # copy and filter the reflections
    reflections = cls._filter_reflections(reflections)

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
      if [beam, goniometer, detector, scan].count(None) != 4:
        raise Sorry('You have a sweep but it is missing a beam, goniometer, '
                    'detector or scan model')
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
      raise Sorry("Goniometer provided, but not enough information about a scan!")
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

    # copy and filter the reflections
    reflections = cls._filter_reflections(reflections)

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
    if nrefs_wo_s1 > 0:
      debug("Setting scattering vectors for %d reflections", nrefs_wo_s1)
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

    # Currently a refinement job can only have one parameterisation of the
    # prediction equation. This can either be of the XYDelPsi (stills) type, the
    # XYPhi (scans) type or the scan-varying XYPhi type with a varying crystal
    # model
    single_as_still = params.refinement.parameterisation.treat_single_image_as_still
    exps_are_stills = []
    for exp in experiments:
      if exp.scan is None:
        exps_are_stills.append(True)
      elif exp.scan.get_num_images() == 1:
        if single_as_still:
          exps_are_stills.append(True)
        elif exp.scan.get_oscillation()[1] == 0.0:
          exps_are_stills.append(True)
        else:
          exps_are_stills.append(False)
      else:
        if exp.scan.get_oscillation()[1] <= 0.0:
          raise Sorry('Cannot refine a zero-width scan')
        exps_are_stills.append(False)

    # check experiment types are consistent
    if not all(exps_are_stills[0] == e for e in exps_are_stills):
      raise Sorry('Cannot refine a mixture of stills and scans')
    do_stills = exps_are_stills[0]

    debug("\nBuilding reflection manager")
    debug("Input reflection list size = %d observations", len(reflections))

    # create reflection manager
    refman = cls.config_refman(params, reflections, experiments, do_stills, verbosity)

    debug("Number of observations that pass initial inclusion criteria = %d",
          refman.get_accepted_refs_size())
    sample_size = refman.get_sample_size()
    if sample_size:
      debug("Working set size = %d observations", sample_size)
    debug("Reflection manager built\n")

    # configure use of sparse data types
    params = cls.config_sparse(params, experiments)

    debug("Building target function")

    # create target function
    target = cls.config_target(params, experiments, refman, do_stills)

    debug("Target function built")

    # create parameterisations
    pred_param, param_reporter, restraints_parameterisation = \
            cls.config_parameterisation(params, experiments, refman, do_stills)

    debug("Prediction equation parameterisation built")
    debug("Parameter order : name mapping")
    for i, e in enumerate(pred_param.get_param_names()):
      debug("Parameter %03d : %s", i, e)

    # Set the prediction equation and restraints parameterisations
    # in the target object
    target.set_prediction_parameterisation(pred_param)
    target.set_restraints_parameterisation(restraints_parameterisation)

    debug("Building refinement engine")

    # create refinery
    refinery = cls.config_refinery(params, target, pred_param, verbosity)

    debug("Refinement engine built")

    # build refiner interface and return
    return Refiner(reflections, experiments, crystal_ids,
                    pred_param, param_reporter, refman, target, refinery,
                    verbosity=verbosity)

  @staticmethod
  def config_sparse(params, experiments):
    """Configure whether to use sparse datatypes"""
    # Automatic selection for sparse parameter
    if params.refinement.parameterisation.sparse == libtbx.Auto:
      if len(experiments) > 1:
        params.refinement.parameterisation.sparse = True
      else:
        params.refinement.parameterisation.sparse = False
      if params.refinement.mp.nproc > 1:
        # sparse vectors cannot be pickled, so can't use easy_mp here
        params.refinement.parameterisation.sparse = False
    # Check incompatible selection
    elif params.refinement.parameterisation.sparse and \
      params.refinement.mp.nproc > 1:
        warning("Could not set sparse=True and nproc={0}".format(
          params.refinement.mp.nproc))
        warning("Resetting sparse=False")
        params.refinement.parameterisation.sparse = False
    return params

  @classmethod
  def config_parameterisation(cls, params, experiments, refman, do_stills):
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
    auto_reduction = params.refinement.parameterisation.auto_reduction

    # Shorten module paths
    import dials.algorithms.refinement.parameterisation as par

    # Get the working set of reflections
    reflections = refman.get_matches()

    # Parameterise unique Beams
    beam_params = []
    for beam in experiments.beams():
      # The Beam is parameterised with reference to a goniometer axis (or None).
      # Use the first (if any) Goniometers this Beam is associated with.
      exp_ids = experiments.indices(beam)
      assoc_gonios = [experiments[i].goniometer for i in exp_ids]
      goniometer = assoc_gonios[0]

      # Parameterise, passing the goniometer (but accepts None)
      beam_param = par.BeamParameterisation(beam, goniometer,
                                                       experiment_ids=exp_ids)
      if beam_options.fix:
        fix_list = [False, False, False]
        if "all" in beam_options.fix:
          fix_list = [True, True, True]
        if "in_spindle_plane" in beam_options.fix:
          fix_list[0] = True
        if "out_spindle_plane" in beam_options.fix:
          fix_list[1] = True
        if "wavelength" in beam_options.fix:
          fix_list[2] = True
        beam_param.set_fixed(fix_list)

      if beam_options.fix_list:
        to_fix = [True if i in beam_options.fix_list else False \
                  for i in range(beam_param.num_total())]
        beam_param.set_fixed(to_fix)

      if beam_param.num_free() > 0:
        beam_params.append(beam_param)

    # Parameterise unique Crystals
    xl_ori_params = []
    xl_uc_params = []
    for crystal in experiments.crystals():
      # This crystal can only ever appear either in scans or in stills
      # (otherwise it requires a different crystal model)
      exp_ids = experiments.indices(crystal)
      assoc_models = [(experiments[i].goniometer, experiments[i].scan) \
                      for i in exp_ids]
      goniometer, scan = assoc_models[0]
      if goniometer is None:
        if not all(g is None and s is None for (g, s) in assoc_models):
          raise Sorry('A crystal model appears in a mixture of scan and still '
                      'experiments, which is not supported')

      if crystal_options.scan_varying:
        # If a crystal is scan-varying, then it must always be found alongside
        # the same Scan and Goniometer in any Experiments in which it appears
        if [goniometer, scan].count(None) != 0:
          raise Sorry('A scan-varying crystal model cannot be created because '
                      'a scan or goniometer model is missing')
        if not all(g is goniometer and s is scan for (g, s) in assoc_models):
          raise Sorry('A single scan-varying crystal model cannot be refined '
                      'when associated with more than one scan or goniometer')

        if crystal_options.num_intervals == "fixed_width":
          sweep_range_deg = scan.get_oscillation_range(deg=True)
          deg_per_interval = crystal_options.interval_width_degrees
          n_intervals = max(int(
            abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)
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

      # get number of fixable units, either parameters or parameter sets in
      # the scan-varying case
      try:
        num_ori = xl_ori_param.num_sets()
      except AttributeError:
        num_ori = xl_ori_param.num_total()
      try:
        num_uc = xl_uc_param.num_sets()
      except AttributeError:
        num_uc = xl_uc_param.num_total()

      if crystal_options.fix:
        if crystal_options.fix == "all":
          xl_ori_param.set_fixed([True] * num_ori)
          xl_uc_param.set_fixed([True] * num_uc)
        elif crystal_options.fix == "cell":
          xl_uc_param.set_fixed([True] * num_uc)
        elif crystal_options.fix == "orientation":
          xl_ori_param.set_fixed([True] * num_ori)
        else: # can only get here if refinement.phil is broken
          raise RuntimeError("crystal_options.fix value not recognised")

      if crystal_options.unit_cell.fix_list:
        to_fix = [True if i in crystal_options.cell.fix_list else False \
                  for i in range(num_uc)]
        xl_uc_param.set_fixed(to_fix)

      if crystal_options.orientation.fix_list:
        to_fix = [True if i in crystal_options.orientation.fix_list else False \
                  for i in range(num_ori)]
        xl_ori_param.set_fixed(to_fix)

      if xl_ori_param.num_free() > 0:
        xl_ori_params.append(xl_ori_param)
      if xl_uc_param.num_free() > 0:
        xl_uc_params.append(xl_uc_param)

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

      if det_param.num_free() > 0:
        det_params.append(det_param)

    # Parameter auto reduction options
    def model_nparam_minus_nref(p, reflections):
      exp_ids = p.get_experiment_ids()
      # Do we have enough reflections to support this parameterisation?
      nparam = p.num_free()
      cutoff = auto_reduction.min_nref_per_parameter * nparam
      isel = flex.size_t()
      for exp_id in exp_ids:
        isel.extend((reflections['id'] == exp_id).iselection())
      nref = len(isel)
      return nref - cutoff

    def panel_gp_nparam_minus_nref(p, pnl_ids, group, reflections):
      exp_ids = p.get_experiment_ids()
      # Do we have enough reflections to support this parameterisation?
      gp_params = [gp == group for gp in p.get_param_panel_groups()]
      fixlist = p.get_fixed()
      free_gp_params = [a and not b for a,b in zip(gp_params, fixlist)]
      nparam = free_gp_params.count(True)
      cutoff = auto_reduction.min_nref_per_parameter * nparam
      isel = flex.size_t()
      for exp_id in exp_ids:
        subsel = (reflections['id'] == exp_id).iselection()
        panels = reflections['panel'].select(subsel)
        for pnl in pnl_ids:
          isel.extend(subsel.select(panels == pnl))
      nref = len(isel)
      return nref - cutoff

    def weak_parameterisation_search(
      beam_params, xl_ori_params, xl_uc_params, det_params, reflections):
      weak = None
      nref_deficit = 0
      panels = None
      pnl_gp = None
      name = None
      for i, p in enumerate(beam_params):
        net_nref = model_nparam_minus_nref(p, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          name = 'Beam{0}'.format(i)
      for i, p in enumerate(xl_ori_params):
        net_nref = model_nparam_minus_nref(p, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          name = 'Crystal{0} orientation'.format(i)
      for i, p in enumerate(xl_uc_params):
        net_nref = model_nparam_minus_nref(p, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          name = 'Crystal{0} unit cell'.format(i)
      for i, p in enumerate(det_params):
        try:
          pnl_groups = p.get_panel_ids_by_group()
          for igp, gp in enumerate(pnl_groups):
            net_nref = panel_gp_nparam_minus_nref(p, gp, igp, reflections)
            if net_nref < nref_deficit:
              nref_deficit = net_nref
              weak = p
              panels = gp
              pnl_gp = igp
              name = 'Detector{0}PanelGroup{1}'.format(i, pnl_gp)
        except:
          net_nref = model_nparam_minus_nref(p, reflections)
          if net_nref < nref_deficit:
            nref_deficit = net_nref
            weak = p
            panels = None
            pnl_gp = None
            name = 'Detector{0}'.format(i)
      return {'parameterisation':weak,
              'panels':panels,
              'panel_group_id':pnl_gp,
              'name':name}

    if auto_reduction.action == 'fail':
      failmsg = 'Too few reflections to parameterise {0}'
      failmsg += '\nTry modifying refinement.parameterisation.auto_reduction options'
      for i, bp in enumerate(beam_params):
        if model_nparam_minus_nref(bp, reflections) < 0:
          mdl = 'Beam{0}'.format(i)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

      for i, xlo in enumerate(xl_ori_params):
        if model_nparam_minus_nref(xlo, reflections) < 0:
          mdl = 'Crystal{0} orientation'.format(i)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

      for i, xluc in enumerate(xl_uc_params):
        if model_nparam_minus_nref(xluc, reflections) < 0:
          mdl = 'Crystal{0} unit cell'.format(i)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

      for i, dp in enumerate(det_params):
        try: # test for hierarchical detector parameterisation
          pnl_groups = dp.get_panel_ids_by_group()
          for igp, gp in enumerate(pnl_groups):
            if panel_gp_nparam_minus_nref(dp, gp, igp, reflections) < 0:
              msg = 'Too few reflections to parameterise Detector{0} panel group {1}'
              msg = msg.format(i, igp)
              msg += '\nTry modifying refinement.parameterisation.auto_reduction options'
              raise Sorry(msg)
        except AttributeError:
          if model_nparam_minus_nref(dp, reflections) < 0:
            mdl = 'Detector{0}'.format(i)
            msg = failmsg.format(mdl)
            raise Sorry(msg)

    elif auto_reduction.action == 'fix':
      warnmsg = 'Too few reflections to parameterise {0}'
      tmp = []
      for i, bp in enumerate(beam_params):
        if model_nparam_minus_nref(bp, reflections) >= 0:
          tmp.append(bp)
        else:
          mdl = 'Beam{0}'.format(i)
          msg = warnmsg.format(mdl)
          warning(msg)
      beam_params = tmp

      tmp = []
      for i, xlo in enumerate(xl_ori_params):
        if model_nparam_minus_nref(xlo, reflections) >= 0:
          tmp.append(xlo)
        else:
          mdl = 'Crystal{0} orientation'.format(i)
          msg = warnmsg.format(mdl)
          warning(msg)
      xl_ori_params = tmp

      tmp = []
      for i, xluc in enumerate(xl_uc_params):
        if model_nparam_minus_nref(xluc, reflections) >= 0:
          tmp.append(xluc)
        else:
          mdl = 'Crystal{0} unit cell'.format(i)
          msg = warnmsg.format(mdl)
          warning(msg)
      xl_uc_params = tmp

      tmp = []
      for i, dp in enumerate(det_params):
        fixlist = dp.get_fixed()
        try: # test for hierarchical detector parameterisation
          pnl_groups = dp.get_panel_ids_by_group()
          for igp, gp in enumerate(pnl_groups):
            if panel_gp_nparam_minus_nref(dp, gp, igp, reflections) < 0:
              msg = 'Too few reflections to parameterise Detector{0}PanelGroup{1}'
              msg = msg.format(i, igp)
              warning(msg)
              gp_params = [gp == igp for gp in dp.get_param_panel_groups()]
              for j, val in enumerate(gp_params):
                if val: fixlist[j] = True
          dp.set_fixed(fixlist)
          if dp.num_free() > 0:
            tmp.append(dp)
          else:
            msg = 'No parameters remain free for Detector{0}'.format(i)
            warning(msg)
        except AttributeError:
          if model_nparam_minus_nref(dp, reflections) >= 0:
            tmp.append(dp)
          else:
            mdl = 'Detector{0}'.format(i)
            msg = warnmsg.format(mdl)
            warning(msg)
      det_params = tmp

    elif auto_reduction.action == 'remove':
      warnmsg = 'Too few reflections to parameterise {0}'
      warnmsg += '\nAssociated reflections will be removed from the Reflection Manager'
      while True:
        dat = weak_parameterisation_search(beam_params,
          xl_ori_params, xl_uc_params, det_params, reflections)
        if dat['parameterisation'] is None: break
        exp_ids = dat['parameterisation'].get_experiment_ids()
        if dat['panels'] is not None:
          msg = warnmsg.format(dat['name'])
          fixlist = dat['parameterisation'].get_fixed()
          pnl_gps = dat['parameterisation'].get_param_panel_groups()
          for i, gp in enumerate(pnl_gps):
            if gp == dat['panel_group_id']: fixlist[i] = True
          dat['parameterisation'].set_fixed(fixlist)
          # identify observations on this panel group from associated experiments
          obs = refman.get_obs()
          isel=flex.size_t()
          for exp_id in exp_ids:
            subsel = (obs['id'] == exp_id).iselection()
            panels_this_exp = obs['panel'].select(subsel)
            for pnl in dat['panels']:
              isel.extend(subsel.select(panels_this_exp == pnl))
        else:
          msg = warnmsg.format(dat['name'])
          fixlist = [True] * dat['parameterisation'].num_total()
          dat['parameterisation'].set_fixed(fixlist)
          # identify observations from the associated experiments
          obs = refman.get_obs()
          isel=flex.size_t()
          for exp_id in exp_ids:
            isel.extend((obs['id'] == exp_id).iselection())
        # Now remove the selected reflections
        sel = flex.bool(len(obs), True)
        sel.set_selected(isel, False)
        refman.filter_obs(sel)
        reflections = refman.get_matches()
        warning(msg)

      # Strip out parameterisations with zero free parameters
      beam_params = [p for p in beam_params if p.num_free() > 0]
      xl_ori_params = [p for p in xl_ori_params if p.num_free() > 0]
      xl_uc_params = [p for p in xl_uc_params if p.num_free() > 0]
      det_params = [p for p in det_params if p.num_free() > 0]

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
        from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import StillsPredictionParameterisationSparse as StillsPredictionParameterisation
      else:
        from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import StillsPredictionParameterisation
      pred_param = StillsPredictionParameterisation(
          experiments,
          det_params, beam_params, xl_ori_params, xl_uc_params)

    else: # doing scans
      if crystal_options.scan_varying:
        if crystal_options.UB_model_per == "reflection":
          if sparse:
            from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
              import VaryingCrystalPredictionParameterisationSparse as PredParam
          else:
            from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
              import VaryingCrystalPredictionParameterisation as PredParam
        elif crystal_options.UB_model_per == "image" or crystal_options.UB_model_per == "block":
          if sparse:
            from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
              import VaryingCrystalPredictionParameterisationFastSparse as PredParam
          else:
            from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
              import VaryingCrystalPredictionParameterisationFast as PredParam
        else:
          raise RuntimeError("UB_model_per=" + crystal_options.scan_varying +
                             " is not a recognised option")
        pred_param = PredParam(
              experiments,
              det_params, beam_params, xl_ori_params, xl_uc_params)
      else:
        if sparse:
          from dials.algorithms.refinement.parameterisation.prediction_parameters \
            import XYPhiPredictionParameterisationSparse as PredParam
        else:
          from dials.algorithms.refinement.parameterisation.prediction_parameters \
            import XYPhiPredictionParameterisation as PredParam
        pred_param = PredParam(
            experiments,
            det_params, beam_params, xl_ori_params, xl_uc_params)

    # Parameter reporting
    param_reporter = par.ParameterReporter(det_params, beam_params,
                                           xl_ori_params, xl_uc_params)

    return pred_param, param_reporter, restraints_param

  @staticmethod
  def config_restraints(params, det_params, beam_params,
        xl_ori_params, xl_uc_params):
    """Given a set of user parameters plus model parameterisations, create a
    restraints plus a parameterisation of these restraints

    Params:
        params The input parameters
        det_params A list of detector parameterisations
        beam_params A list of beam parameterisations,
        xl_ori_params A list of crystal orientation parameterisations
        xl_uc_params A list of crystal unit cell parameterisations

    Returns:
        A restraints parameterisation
    """

    from dials.algorithms.refinement.restraints import RestraintsParameterisation
    rp = RestraintsParameterisation(detector_parameterisations = det_params,
               beam_parameterisations = beam_params,
               xl_orientation_parameterisations = xl_ori_params,
               xl_unit_cell_parameterisations = xl_uc_params)

    # Shorten params path
    # FIXME Only unit cell restraints currently supported
    #beam_r = params.refinement.parameterisation.beam.restraints
    cell_r = params.refinement.parameterisation.crystal.unit_cell.restraints
    #orientation_r = params.refinement.parameterisation.crystal.orientation.restraints
    #detector_r = params.refinement.parameterisation.detector.restraints

    for tie in cell_r.tie_to_target:
      if len(tie.values) != 6:
        raise Sorry("6 cell parameters must be provided as the tie_to_target.values.")
      if len(tie.sigmas) != 6:
        raise Sorry("6 sigmas must be provided as the tie_to_target.sigmas. "
                    "Note that individual sigmas of 0.0 will remove "
                    "the restraint for the corresponding cell parameter.")
      if not tie.id:
        raise Sorry("At least one experiment id must be provided as the "
                    "tie_to_target.id")
      for exp_id in tie.id:
        rp.add_restraints_to_target_xl_unit_cell(exp_id, tie.values, tie.sigmas)

    # FIXME Group ties not available yet
    for tie in cell_r.tie_to_group:
      pass

    return rp

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

    import dials.algorithms.refinement.engine as engine # implicit import

    if options.engine == "SimpleLBFGS":
      from engine import SimpleLBFGS as refinery
    elif options.engine == "LBFGScurvs":
      from engine import LBFGScurvs as refinery
    elif options.engine == "GaussNewton":
      from engine import GaussNewtonIterations as refinery
    elif options.engine == "LevMar":
      from engine import LevenbergMarquardtIterations as refinery
    elif options.engine == "SparseLevMar":
      from sparse_engine import SparseLevenbergMarquardtIterations as refinery
    else:
      raise RuntimeError("Refinement engine " + options.engine +
                         " not recognised")

    debug("Selected refinement engine type: %s", options.engine)

    engine = refinery(target = target,
            prediction_parameterisation = pred_param,
            log = options.log,
            verbosity = verbosity,
            track_step = options.track_step,
            track_gradient = options.track_gradient,
            track_parameter_correlation = options.track_parameter_correlation,
            track_out_of_sample_rmsd = options.track_out_of_sample_rmsd,
            max_iterations = options.max_iterations)

    if params.refinement.mp.nproc > 1:
      nproc = params.refinement.mp.nproc
      try:
        engine.set_nproc(nproc)
      except NotImplementedError:
        warning("Could not set nproc={0} for refinement engine of type {1}".format(
          nproc, options.engine))

    return engine

  @staticmethod
  def config_refman(params, reflections, experiments, do_stills, verbosity):
    """Given a set of parameters and models, build a reflection manager

    Params:
        params The input parameters

    Returns:
        The reflection manager instance
    """

    # Shorten parameter path
    options = params.refinement.reflections
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
      debug("Random seed set to %d", options.random_seed)

    # check whether we deal with stills or scans
    if do_stills:
      from dials.algorithms.refinement.reflection_manager import \
          StillsReflectionManager as refman
      # check incompatible weighting strategy
      if options.weighting_strategy.override == "statistical":
        raise Sorry('The "statistical" weighting strategy is not compatible '
                    'with stills refinement')
    else:
      from dials.algorithms.refinement.reflection_manager import ReflectionManager as refman
      # check incompatible weighting strategy
      if options.weighting_strategy.override == "stills":
        raise Sorry('The "stills" weighting strategy is not compatible with '
                    'scan refinement')

    # set automatic outlier rejection options
    if options.outlier.algorithm in ('auto', libtbx.Auto):
      if do_stills:
        options.outlier.algorithm = 'sauter_poon'
      else:
        options.outlier.algorithm = 'mcd'

    if options.outlier.separate_panels is libtbx.Auto:
      if do_stills:
        options.outlier.separate_panels = False
      else:
        options.outlier.separate_panels = True

    if options.outlier.algorithm == 'sauter_poon':
      if options.outlier.sauter_poon.px_sz is libtbx.Auto:
        # get this from the first panel of the first detector
        options.outlier.sauter_poon.px_sz = experiments.detectors()[0][0].get_pixel_size()

    # do outlier rejection?
    if options.outlier.algorithm in ("null", None):
      outlier_detector = None
    else:
      if do_stills:
        colnames = ["x_resid", "y_resid"]
      else:
        colnames = ["x_resid", "y_resid", "phi_resid"]
      from dials.algorithms.refinement.outlier_detection import CentroidOutlierFactory
      outlier_detector = CentroidOutlierFactory.from_parameters_and_colnames(
        options, colnames)

    # override default weighting strategy?
    weighting_strategy = None
    if options.weighting_strategy.override == "statistical":
      from dials.algorithms.refinement.weighting_strategies \
        import StatisticalWeightingStrategy
      weighting_strategy = StatisticalWeightingStrategy()
    elif options.weighting_strategy.override == "stills":
      from dials.algorithms.refinement.weighting_strategies \
        import StillsWeightingStrategy
      weighting_strategy = StillsWeightingStrategy(
        options.weighting_strategy.delpsi_constant)
    elif options.weighting_strategy.override == "constant":
      from dials.algorithms.refinement.weighting_strategies \
        import ConstantWeightingStrategy
      weighting_strategy = ConstantWeightingStrategy(
        *options.weighting_strategy.constants, stills=do_stills)

    # calculate reflection block_width?
    if params.refinement.parameterisation.crystal.scan_varying:
      from dials.algorithms.refinement.reflection_manager import BlockCalculator
      block_calculator = BlockCalculator(experiments, reflections)
      if params.refinement.parameterisation.crystal.UB_model_per == "block":
        reflections = block_calculator.per_width(options.block_width, deg=True)
      elif params.refinement.parameterisation.crystal.UB_model_per == "image":
        reflections = block_calculator.per_image()

    return refman(reflections=reflections,
            experiments=experiments,
            nref_per_degree=nref_per_degree,
            max_sample_size = options.maximum_sample_size,
            min_sample_size = options.minimum_sample_size,
            close_to_spindle_cutoff=options.close_to_spindle_cutoff,
            outlier_detector=outlier_detector,
            weighting_strategy_override=weighting_strategy,
            verbosity=verbosity)

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

    # build managed reflection predictors
    from dials.algorithms.refinement.prediction import ExperimentsPredictor
    ref_predictor = ExperimentsPredictor(experiments, do_stills)

    # Determine whether the target is in X, Y, Phi space or just X, Y.
    if do_stills:
      if sparse:
        from dials.algorithms.refinement.target_stills \
          import LeastSquaresStillsResidualWithRmsdCutoffSparse as targ
      else:
        from dials.algorithms.refinement.target_stills \
          import LeastSquaresStillsResidualWithRmsdCutoff as targ
    else:
      if sparse:
        from dials.algorithms.refinement.target \
          import LeastSquaresPositionalResidualWithRmsdCutoffSparse as targ
      else:
        from dials.algorithms.refinement.target \
          import LeastSquaresPositionalResidualWithRmsdCutoff as targ

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
    get_matches
    get_param_reporter
    parameter_correlation_plot
    selection_used_for_refinement
    predict_for_reflection_table
    predict_for_indexed

  Notes:
    * The return value of run is a recorded history of the refinement
    * The experiments accessor provides a copy of the experiments used by
      refinement
    * The model accessors provide copies of those models that might be modified
      by refinement (beam, crystal and detector) TO BE DEPRECATED
    * get_matches exposes the function of the same name from the privately
      stored reflection manager
    * The return value of selection_used_for_refinement is a flex.bool

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

  def rmsds(self):
    """Return rmsds of the current model"""

    # ensure predictions for the matches are up to date
    self._refinery.prepare_for_step()

    return self._target.rmsds()

  def rmsds_for_reflection_table(self, reflections):
    """Calculate unweighted RMSDs for the specified reflections"""

    # ensure predictions for these reflections are up to date
    preds = self.predict_for_reflection_table(reflections)

    return self._target.rmsds_for_reflection_table(preds)


  def get_matches(self):
    """Delegated to the reflection manager"""

    # FIXME Consider: Does this information really need to be exposed by the
    # public API (indexing code seems to use it, but is it necessary?)
    return self._refman.get_matches()

  def get_free_reflections(self):
    """Delegated to the reflection manager"""

    return self._refman.get_free_reflections()

  def get_param_reporter(self):
    """Get the ParameterReport object linked to this Refiner"""

    return self._param_report

  def get_parameter_correlation_matrix(self, step, col_select=None):
    """Return the correlation matrix between columns of the Jacobian at
    the specified refinement step. The parameter col_select can be used
    to select subsets of the full number of columns. The column labels
    are also returned as a list of strings"""

    corrmat = self._refinery.get_correlation_matrix_for_step(step)
    if corrmat is None: return None
    assert corrmat.is_square_matrix()

    if col_select is None:
      col_select = range(corrmat.all()[0])

    all_labels = self._pred_param.get_param_names()
    idx = []
    for col in col_select:
      try: # column is specified by name
        idx.append(all_labels.index(col))
      except ValueError: # column specified by number
        try:
          idx.append(int(col))
        except ValueError:
          info("Invalid selection of columns for correlation plot. " + \
               "No plot will be produced")
          return None
    labels = [all_labels[e] for e in idx]
    num_cols = num_rows = len(labels)

    from scitbx.array_family import flex
    sub_corrmat = flex.double(flex.grid(num_cols, num_cols))

    for (i, x) in enumerate(idx):
      for (j, y) in enumerate(idx):
        sub_corrmat[i,j] = corrmat[x, y]

    return (sub_corrmat, labels)

  def print_step_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    info("\nRefinement steps:")

    rmsd_multipliers = []
    header = ["Step", "Nref"]
    for (name, units) in zip(self._target.rmsd_names, self._target.rmsd_units):
      if units == "mm":
        header.append(name + "\n(mm)")
        rmsd_multipliers.append(1.0)
      elif units == "rad": # convert radians to degrees for reporting
        header.append(name + "\n(deg)")
        rmsd_multipliers.append(rad2deg)
      else: # leave unknown units alone
        header.append(name + "\n(" + units + ")")

    rows = []
    for i in range(self._refinery.history.get_nrows()):
      rmsds = [r*m for (r,m) in zip(self._refinery.history["rmsd"][i], rmsd_multipliers)]
      rows.append([str(i), str(self._refinery.history["num_reflections"][i])] + \
        ["%.5g" % r for r in rmsds])

    st = simple_table(rows, header)
    info(st.format())
    info(self._refinery.history.reason_for_termination)

    return

  def print_out_of_sample_rmsd_table(self):
    """print out-of-sample RSMDs per step, if these were tracked"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    # check if it makes sense to proceed
    if not self._refinery.history.has_key("out_of_sample_rmsd"): return
    nref = len(self.get_free_reflections())
    if nref < 10: return # don't do anything if very few refs

    info("\nRMSDs for out-of-sample (free) reflections:")

    rmsd_multipliers = []
    header = ["Step", "Nref"]
    for (name, units) in zip(self._target.rmsd_names, self._target.rmsd_units):
      if units == "mm":
        header.append(name + "\n(mm)")
        rmsd_multipliers.append(1.0)
      elif units == "rad": # convert radians to degrees for reporting
        header.append(name + "\n(deg)")
        rmsd_multipliers.append(rad2deg)
      else: # leave unknown units alone
        header.append(name + "\n(" + units + ")")

    rows = []
    for i in range(self._refinery.history.get_nrows()):
      rmsds = [r*m for r, m in zip(self._refinery.history["out_of_sample_rmsd"][i],
                                   rmsd_multipliers)]
      rows.append([str(i), str(nref)] + ["%.5g" % e for e in rmsds])

    st = simple_table(rows, header)
    info(st.format())

    return

  def print_exp_rmsd_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    info("\nRMSDs by experiment:")

    header = ["Exp", "Nref"]
    for (name, units) in zip(self._target.rmsd_names, self._target.rmsd_units):
      if name == "RMSD_X" or name == "RMSD_Y" and units == "mm":
        header.append(name + "\n(px)")
      elif name == "RMSD_Phi" and units == "rad": # convert radians to images for reporting of scans
        header.append("RMSD_Z" + "\n(images)")
      elif name == "RMSD_DeltaPsi" and units == "rad": # convert radians to degrees for reporting of stills
        header.append(name + "\n(deg)")
      else: # skip RMSDs that cannot be expressed in image/scan space
        pass

    rows = []
    for iexp, exp in enumerate(self._experiments):
      detector = exp.detector
      px_sizes = [p.get_pixel_size() for p in detector]
      it = iter(px_sizes)
      px_size = next(it)
      if not all(tst == px_size for tst in it):
        info("The detector in experiment %d does not have the same pixel " + \
             "sizes on each panel. Skipping...", iexp)
        continue
      px_per_mm = [1./e for e in px_size]

      scan = exp.scan
      try:
        temp = scan.get_oscillation(deg=False)
        images_per_rad  = 1./abs(scan.get_oscillation(deg=False)[1])
      except AttributeError:
        images_per_rad = None

      raw_rmsds = self._target.rmsds_for_experiment(iexp)
      if raw_rmsds is None: continue # skip experiments where rmsd cannot be calculated
      num = self._target.get_num_matches_for_experiment(iexp)
      rmsds = []
      for (name, units, rmsd) in zip(self._target.rmsd_names, self._target.rmsd_units, raw_rmsds):
        if name == "RMSD_X" and units == "mm":
          rmsds.append(rmsd * px_per_mm[0])
        elif name == "RMSD_Y" and units == "mm":
          rmsds.append(rmsd * px_per_mm[1])
        elif name == "RMSD_Phi" and units == "rad":
          rmsds.append(rmsd * images_per_rad)
        elif name == "RMSD_DeltaPsi" and units == "rad":
          rmsds.append(rmsd * rad2deg)
      rows.append([str(iexp), str(num)] + ["%.5g" % r for r in rmsds])

    if len(rows) > 0:
      truncated = False
      max_rows = 20
      if self._verbosity < 2 and len(rows) > max_rows:
        rows = rows[0:max_rows]
        truncated = True
      st = simple_table(rows, header)
      info(st.format())
      if truncated:
        info("Table truncated to show the first %d experiments only", max_rows)
        info("Re-run with verbosity >= 2 to show all experiments")

    return

  def print_panel_rmsd_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    info("\nRMSDs by panel:")

    header = ["Panel", "Nref"]
    for (name, units) in zip(self._target.rmsd_names, self._target.rmsd_units):
      if name == "RMSD_X" or name == "RMSD_Y" and units == "mm":
        header.append(name + "\n(px)")
      elif name == "RMSD_Phi" and units == "rad": # convert radians to images for reporting of scans
        header.append("RMSD_Z" + "\n(images)")
      elif name == "RMSD_DeltaPsi" and units == "rad": # convert radians to degrees for reporting of stills
        header.append(name + "\n(deg)")
      else: # skip RMSDs that cannot be expressed in image/scan space
        pass

    rows = []
    for ipanel, panel in enumerate(self._detector):

      px_size = panel.get_pixel_size()
      px_per_mm = [1./e for e in px_size]

      scan = self._scan
      try:
        temp = scan.get_oscillation(deg=False)
        images_per_rad  = 1./abs(scan.get_oscillation(deg=False)[1])
      except AttributeError:
        images_per_rad = None

      num = self._target.get_num_matches_for_panel(ipanel)
      if num <= 0: continue
      raw_rmsds = self._target.rmsds_for_panel(ipanel)
      if raw_rmsds is None: continue # skip panels where rmsd cannot be calculated
      rmsds = []
      for (name, units, rmsd) in zip(self._target.rmsd_names, self._target.rmsd_units, raw_rmsds):
        if name == "RMSD_X" and units == "mm":
          rmsds.append(rmsd * px_per_mm[0])
        elif name == "RMSD_Y" and units == "mm":
          rmsds.append(rmsd * px_per_mm[1])
        elif name == "RMSD_Phi" and units == "rad":
          rmsds.append(rmsd * images_per_rad)
        elif name == "RMSD_DeltaPsi" and units == "rad":
          rmsds.append(rmsd * rad2deg)
      rows.append([str(ipanel), str(num)] + ["%.5g" % r for r in rmsds])

    if len(rows) > 0:
      st = simple_table(rows, header)
      info(st.format())

    return

  def run(self):
    """Run refinement"""

    ####################################
    # Do refinement and return history #
    ####################################

    debug("\nExperimental models before refinement:")
    debug(str(self._beam))
    debug(str(self._detector))
    if self._goniometer: debug(str(self._goniometer))
    if self._scan: debug(str(self._scan))
    for i, x in zip(self._crystal_ids, self._crystals):
      msg = "%d " % i
      msg += str(x)
      debug(msg)

    self._refinery.run()

    self.print_step_table()
    self.print_out_of_sample_rmsd_table()
    self.print_exp_rmsd_table()

    if len(self._detector) > 1:
      self.print_panel_rmsd_table()

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

        # get state covariance matrices the whole range of images. We select
        # the first element of this at each image because crystal scan-varying
        # parameterisations are not multi-state
        state_cov_list = [self._pred_param.calculate_model_state_uncertainties(
          obs_image_number=t, experiment_id=iexp) for t in range(ar_range[0],
                                                            ar_range[1]+1)]
        u_cov_list, b_cov_list = zip(*state_cov_list)

        # return these to the model parameterisations to be set in the models
        self._pred_param.set_model_state_uncertainties(
          u_cov_list, b_cov_list, iexp)

    debug("\nExperimental models after refinement:")
    debug(str(self._beam))
    debug(str(self._detector))
    if self._goniometer: debug(str(self._goniometer))
    if self._scan: debug(str(self._scan))
    for i, x in zip(self._crystal_ids, self._crystals):
      msg = "%d " % i
      msg += str(x)
      debug(msg)

    # Report on the refined parameters
    debug(str(self._param_report))

    # Return the refinement history
    return self._refinery.history

  def selection_used_for_refinement(self):
    """Return a selection as a flex.bool in terms of the input reflection
    data of those reflections that were used in the final step of
    refinement."""

    from scitbx.array_family import flex
    matches = self._refman.get_matches()
    selection = flex.bool(len(self._refman.get_indexed()), False)

    try: # new reflection table format for matches
      isel = matches['iobs']
      selection.set_selected(isel, True)
    except TypeError: # old ObsPredMatch format for matches
      for m in matches:
        selection[m.iobs] = True

    return selection

  def predict_for_indexed(self):
    """perform prediction for all the indexed reflections passed into
    refinement"""

    return self.predict_for_reflection_table(self._refman.get_indexed())

  def predict_for_reflection_table(self, reflections):
    """perform prediction for all reflections in the supplied table"""

    # delegate to the target object, which has access to the predictor
    return self._target.predict_for_reflection_table(reflections)

