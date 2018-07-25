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

from __future__ import absolute_import, division
import logging
logger = logging.getLogger(__name__)

from dxtbx.model.experiment_list import ExperimentList
from dials.array_family import flex
from dials.algorithms.refinement.refinement_helpers import ordinal_number
from libtbx.phil import parse
from libtbx.utils import Sorry
import libtbx

from dials_refinement_helpers_ext import pgnmn_iter as pgnmn
from dials_refinement_helpers_ext import ucnmn_iter as ucnmn
from dials_refinement_helpers_ext import mnmn_iter as mnmn

# The include scope directive does not work here. For example:
#
#   include scope dials.algorithms.refinement.outlier_detection.phil_scope
#
# results in:
#
#   AttributeError: 'module' object has no attribute 'refinement'
#
# to work around this, just include external phil scopes as strings
from dials.algorithms.refinement.outlier_detection.outlier_base \
  import phil_str as outlier_phil_str
from dials.algorithms.refinement.restraints.restraints_parameterisation \
  import uc_phil_str as uc_restraints_phil_str
from dials.algorithms.refinement.constraints import phil_str as constr_phil_str
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
  import phil_str as sv_phil_str
from dials.algorithms.refinement.engine import refinery_phil_str
format_data = {'outlier_phil':outlier_phil_str,
               'uc_restraints_phil':uc_restraints_phil_str,
               'constr_phil':constr_phil_str,
               'sv_phil_str':sv_phil_str,
               'refinery_phil':refinery_phil_str}
phil_scope = parse('''

refinement
  .help = "Parameters to configure the refinement"
{

  mp
    .expert_level = 2
  {
    nproc = 1
      .type = int(value_min=1)
      .help = "The number of processes to use. Not all choices of refinement"
              "engine support nproc > 1. Where multiprocessing is possible,"
              "it is helpful only in certain circumstances, so this is not"
              "recommended for typical use."
  }

  verbosity = 0
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
                "remove these reflections from other parameterisations of the"
                "global model too. For example, if a crystal model could not"
                "be parameterised it will be excised completely and not"
                "contribute to the joint refinement of the detector and beam."
                "In the fix mode, reflections emanating from that crystal will"
                "still form residuals and will contribute to detector and beam"
                "refinement."
        .type = choice

      detector_reduce = False
        .type = bool
        .help = "Special case designed for detector metrology refinement"
                "(particularly of the CSPAD). See detector_reduce_list for"
                "details."
        .expert_level = 1

      detector_reduce_list = Dist Tau2 Tau3
        .type = strings
        .help = "Partial names to match to detector parameters to try fixing."
                "If there are still not"
                "enough parameters for refinement after fixing these, then"
                "fail. This is to ensure that metrology refinement never"
                "completes if it is not able to refine some panels. The default"
                "is to try fixing the distance as well as Tau2 and Tau3"
                "rotations of detector panel, leaving the in-plane shifts and"
                "the rotation around the detector normal for refinement."
                "groups only."
        .expert_level = 1
    }

    scan_varying = False
      .help = "Allow models that are not forced to be static to vary during the
               scan"
      .type = bool

    compose_model_per = image *block
      .help = "For scan-varying parameterisations, compose a new model either"
              "every image or within blocks of a width specified in the"
              "reflections parameters. When this block width is larger than the"
              "image width the result is faster, with a trade-off in accuracy"
      .type = choice
      .expert_level = 1

    debug_centroid_analysis = False
      .help = "Set True to write out a file containing the reflections used"
              "for centroid analysis for automatic setting of the  scan-varying"
              "interval width. This can then be analysed with"
              "dev.dials.plot_centroid_analysis"
      .type = bool
      .expert_level = 2

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
        .type = strings
        .help = "Fix specified parameters by a list of 0-based indices or"
                "partial names to match"
        .expert_level = 1

      %(constr_phil)s

      force_static = True
        .type = bool
        .help = "Force a static parameterisation for the beam when doing"
                "scan-varying refinement"
        .expert_level = 1

      %(sv_phil_str)s
    }

    crystal
      .help = "crystal parameters"
    {
      fix = all cell orientation
        .help = "Fix crystal parameters"
        .type = choice

      unit_cell
        .expert_level = 1
      {
        fix_list = None
          .type = strings
          .help = "Fix specified parameters by a list of 0-based indices or"
                  "partial names to match"
          .expert_level = 1

        %(uc_restraints_phil)s

        %(constr_phil)s

        force_static = False
          .type = bool
          .help = "Force a static parameterisation for the crystal unit cell"
                  "when doing scan-varying refinement"
          .expert_level = 1

        set_scan_varying_errors = False
          .type = bool
          .help = "If scan-varying refinement is done, and if the estimated"
                  "covariance of the B matrix has been calculated by the"
                  "minimiser, choose whether to return this to the model or"
                  "not. The default is not to, in order to keep the file size"
                  "of the serialized model small."

        %(sv_phil_str)s
      }

      orientation
        .expert_level = 1
      {
        fix_list = None
          .type = strings
          .help = "Fix specified parameters by a list of 0-based indices or"
                  "partial names to match"
          .expert_level = 1

        %(constr_phil)s

        force_static = False
          .type = bool
          .help = "Force a static parameterisation for the crystal orientation"
                  "when doing scan-varying refinement"
          .expert_level = 1

        %(sv_phil_str)s
      }
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
        .type = strings
        .help = "Fix specified parameters by a list of 0-based indices or"
                "partial names to match"
        .expert_level = 1

      %(constr_phil)s

      force_static = True
        .type = bool
        .help = "Force a static parameterisation for the detector"
                "when doing scan-varying refinement"
        .expert_level = 1

      %(sv_phil_str)s
    }

    goniometer
      .help = "goniometer setting matrix parameters"
    {
      fix = *all in_beam_plane out_beam_plane
        .help = "Whether to fix goniometer parameters. By default,"
                "fix all. Alternatively the setting matrix can be constrained"
                "to allow rotation only within the spindle-beam plane"
                "or to allow rotation only around an axis that lies in that"
                "plane. Set to None to refine the in two orthogonal directions."
        .type = choice(multi=True)

      fix_list = None
        .type = strings
        .help = "Fix specified parameters by a list of 0-based indices or"
                "partial names to match"
        .expert_level = 1

      %(constr_phil)s

      force_static = True
        .type = bool
        .help = "Force a static parameterisation for the goniometer when doing"
                "scan-varying refinement"
        .expert_level = 1

      %(sv_phil_str)s
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

    spherical_relp_model = False
      .help = "For stills refinement, set true to use the spherical relp model"
              "for prediction and gradients."
      .type = bool
      .expert_level = 1
  }

  %(refinery_phil)s

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

    bin_size_fraction = 0.0
      .help = "Set this to a fractional value, say 0.2, to make a cut off in"
              "the natural discrete units of positional data, viz.,"
              "(pixel width, pixel height, image thickness in phi). This would"
              "then determine when the RMSD target is achieved. Only used if"
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

    reflections_per_degree = None
      .help = "The number of centroids per degree of the sweep to use in"
              "refinement. Set to None to use all suitable reflections."
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

    trim_scan_edges = 0.0
      .help = "Reflections within this value in degrees from the centre of the"
              "first or last image of the scan will be removed before"
              "refinement, unless doing so would result in too few remaining"
              "reflections. Reflections that are truncated at the scan edges"
              "have poorly-determined centroids and can bias the refined model"
              "if they are included."
      .type = float(value_min=0,value_max=1)
      .expert_level = 1

    block_width = 1.0
      .help = "Width of a reflection 'block' (in degrees) determining how fine-"
              "grained the model used for scan-varying prediction during"
              "refinement is. Currently only has any effect if the crystal"
              "parameterisation is set to use compose_model_per=block"
      .type = float(value_min = 0.)
      .expert_level = 1

    weighting_strategy
      .help = "Parameters to configure weighting strategy overrides"
      .expert_level = 1
    {
      override = statistical stills constant external_deltapsi
        .help = "selection of a strategy to override default weighting behaviour"
        .type = choice

      delpsi_constant = 1000000
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
            'xyzobs.mm.variance', 'flags', 'delpsical.weights']
    # NB xyzobs.px.value & xyzcal.px required by SauterPoon outlier rejector
    # NB delpsical.weights is used by ExternalDelPsiWeightingStrategy
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
                                       verbosity=None,
                                       copy_experiments=True):

    #TODO Checks on the input
    #E.g. does every experiment contain at least one overlapping model with at
    #least one other experiment? Are all the experiments either rotation series
    #or stills (the combination of both not yet supported)?

    # if no verbosity override is given, take from the parameters
    if verbosity is None:
      verbosity = params.refinement.verbosity

    if copy_experiments:
      # copy the experiments
      import copy
      experiments = copy.deepcopy(experiments)

    # copy and filter the reflections
    reflections = cls._filter_reflections(reflections)

    return cls._build_components(params,
                                 reflections,
                                 experiments,
                                 verbosity=verbosity,
                                 copy_experiments=copy_experiments)

  @classmethod
  def _build_components(cls, params, reflections, experiments,
                        verbosity, copy_experiments=True):
    """low level build"""

    if verbosity == 0:
      logger.disabled = True

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

    logger.debug("\nBuilding reflection manager")
    logger.debug("Input reflection list size = %d observations", len(reflections))

    # create reflection manager
    refman = cls.config_refman(params, reflections, experiments, do_stills, verbosity)

    logger.debug("Number of observations that pass initial inclusion criteria = %d",
          refman.get_accepted_refs_size())
    sample_size = refman.get_sample_size()
    if sample_size > 0:
      logger.debug("Working set size = %d observations", sample_size)
    logger.debug("Reflection manager built\n")

    # configure use of sparse data types
    params = cls.config_sparse(params, experiments)

    logger.debug("Building target function")

    # create target function
    target = cls.config_target(params, experiments, refman, do_stills)

    logger.debug("Target function built")

    # determine whether to do basic centroid analysis to automatically
    # determine outlier rejection block
    if params.refinement.reflections.outlier.block_width is libtbx.Auto:
      ca = refman.get_centroid_analyser()
      analysis = ca(calc_average_residuals=False, calc_periodograms=False)
    else:
      analysis = None

    # Now predictions and centroid analysis are available, so we can finalise
    # the reflection manager
    refman.finalise(analysis)

    # create parameterisations
    pred_param, param_reporter, restraints_parameterisation = \
      cls.config_parameterisation(params, experiments, refman, do_stills)

    logger.debug("Prediction equation parameterisation built")
    logger.debug("Parameter order : name mapping")
    for i, e in enumerate(pred_param.get_param_names()):
      logger.debug("Parameter %03d : %s", i + 1, e)

    # Set the prediction equation and restraints parameterisations
    # in the target object
    target.set_prediction_parameterisation(pred_param)
    target.set_restraints_parameterisation(restraints_parameterisation)

    # Build a constraints manager, if requested
    from dials.algorithms.refinement.constraints import ConstraintManagerFactory
    cmf = ConstraintManagerFactory(params, pred_param, verbosity)
    constraints_manager = cmf()

    logger.debug("Building refinement engine")

    # create refinery
    refinery = cls.config_refinery(params, target, pred_param,
      constraints_manager, verbosity)

    logger.debug("Refinement engine built")

    # build refiner interface and return
    return Refiner(reflections, experiments,
                    pred_param, param_reporter, refman, target, refinery,
                    verbosity=verbosity, copy_experiments=copy_experiments)

  @staticmethod
  def config_sparse(params, experiments):
    """Configure whether to use sparse datatypes"""
    # Automatic selection for sparse parameter
    if params.refinement.parameterisation.sparse == libtbx.Auto:
      if len(experiments) > 1:
        params.refinement.parameterisation.sparse = True
      else:
        params.refinement.parameterisation.sparse = False
      if params.refinement.refinery.engine == "SparseLevMar":
        params.refinement.parameterisation.sparse = True
      if params.refinement.mp.nproc > 1:
        if params.refinement.refinery.engine != "SparseLevMar":
          # sparse vectors cannot be pickled, so can't use easy_mp here
          params.refinement.parameterisation.sparse = False
        else:
          pass # but SparseLevMar requires sparse jacobian; does not implement mp
    # Check incompatible selection
    elif params.refinement.parameterisation.sparse and \
      params.refinement.mp.nproc > 1:
        logger.warning("Could not set sparse=True and nproc={0}".format(
          params.refinement.mp.nproc))
        logger.warning("Resetting sparse=False")
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

    # Shorten options path
    options = params.refinement.parameterisation

    # Shorten module paths
    import dials.algorithms.refinement.parameterisation as par

    # function to convert fix_lists into to_fix selections
    from dials.algorithms.refinement.refinement_helpers import string_sel

    # Get the working set of reflections
    reflections = refman.get_matches()

    # Define a helper function for parameter fixing
    def filter_parameter_names(parameterisation):
      # scan-varying suffixes like '_sample1' should be removed from
      # the list of parameter names so that it is num_free in length
      import re
      pattern = re.compile(r"_sample[0-9]+$")
      names = [pattern.sub('', e) for e in parameterisation.get_param_names(only_free=False)]
      filtered_names = []
      for name in names:
        if name not in filtered_names: filtered_names.append(name)
      return filtered_names

    # If required, do full centroid analysis (now on the outlier-rejected
    # reflections) to determine suitable interval widths for scan-varying
    # refinement
    analysis = None
    if options.scan_varying:
      tst = [options.beam.smoother,
             options.crystal.orientation.smoother,
             options.crystal.unit_cell.smoother,
             options.detector.smoother,
             options.goniometer.smoother]
      tst = [(e.absolute_num_intervals is None and
              e.interval_width_degrees is libtbx.Auto) for e in tst]
      if any(tst):
        logger.info('Doing centroid analysis to '
          'automatically determine scan-varying interval widths')
        ca = refman.get_centroid_analyser(debug=options.debug_centroid_analysis)
        analysis = ca()

    # Use the results of centroid analysis to suggest suitable interval widths
    # for each experiment. This will be the smallest of the proposed intervals
    # for each of the residuals in x, y and phi, as long as this is not smaller
    # than either the outlier rejection block width, or 9.0 degrees.
    if analysis is not None:
      for i, a in enumerate(analysis):
        intervals = [a.get('x_interval'),
                     a.get('y_interval'),
                     a.get('phi_interval')]
        try:
          min_interval = min([e for e in intervals if e is not None])
        except ValueError:
          # empty list - analysis was unable to suggest a suitable interval
          # width. Default to the safest case
          phi_min, phi_max  = experiments[i].scan.get_oscillation_range(deg=True)
          a['interval_width'] = abs(phi_max - phi_min)
          logger.info('Exp id {0} suggested interval width could not be '
              'determined and will be reset to the scan width of '
              '{1:.1f} degrees'.format(i, a['interval_width']))
          continue
        min_interval = max(min_interval, 9.0)
        block_size = a.get('block_size')
        if block_size is not None:
          min_interval = max(min_interval, block_size)
        a['interval_width'] = min_interval
        logger.info('Exp id {0} suggested interval width = {1:.1f} degrees'.format(
            i, min_interval))

    # Parameterise unique Beams
    beam_params = []
    for ibeam, beam in enumerate(experiments.beams()):
      # The Beam is parameterised with reference to a goniometer axis (or None).
      # Use the first (if any) Goniometers this Beam is associated with.
      exp_ids = experiments.indices(beam)
      assoc_models = [(experiments[i].goniometer, experiments[i].scan) \
                      for i in exp_ids]
      goniometer, scan = assoc_models[0]

      if options.scan_varying:
        if not options.beam.force_static:
          # If a beam is scan-varying, then it must always be found alongside
          # the same Scan and Goniometer in any Experiments in which it appears
          if [goniometer, scan].count(None) != 0:
            raise Sorry('A scan-varying beam model cannot be created because '
                        'a scan or goniometer model is missing')
          if not all(g is goniometer and s is scan for (g, s) in assoc_models):
            raise Sorry('A single scan-varying beam model cannot be refined '
                        'when associated with more than one scan or goniometer')
        sweep_range_deg = scan.get_oscillation_range(deg=True)
        array_range = scan.get_array_range()

        if not options.beam.force_static:
          n_intervals = options.beam.smoother.absolute_num_intervals
          if n_intervals is None:
            deg_per_interval = options.beam.smoother.interval_width_degrees
            if deg_per_interval is libtbx.Auto and analysis is not None:
              intervals = [analysis[i]['interval_width'] for i in exp_ids]
              deg_per_interval = min(intervals)
            elif deg_per_interval is None:
              deg_per_interval = 36.0
            n_intervals = max(int(
              abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)

          beam_param = par.ScanVaryingBeamParameterisation(
                                              beam,
                                              array_range,
                                              n_intervals,
                                              goniometer=goniometer,
                                              experiment_ids=exp_ids)
        else: # force model to be static
          beam_param = par.BeamParameterisation(beam, goniometer,
                                                       experiment_ids=exp_ids)
      else:
        # Parameterise scan static beam, passing the goniometer
        beam_param = par.BeamParameterisation(beam, goniometer,
                                                       experiment_ids=exp_ids)

      # get number of fixable units, either parameters or parameter sets in
      # the scan-varying case
      try:
        num_beam = beam_param.num_sets()
      except AttributeError:
        num_beam = beam_param.num_total()

      fix_list = []
      if options.beam.fix_list:
        fix_list.extend(options.beam.fix_list)

      if options.beam.fix:
        if "all" in options.beam.fix:
          beam_param.set_fixed([True] * num_beam)
        if "in_spindle_plane" in options.beam.fix:
          fix_list.append('Mu1')
        if "out_spindle_plane" in options.beam.fix:
          fix_list.append('Mu2')
        if "wavelength" in options.beam.fix:
          fix_list.append('nu')

      if fix_list:
        names = filter_parameter_names(beam_param)
        assert len(names) == num_beam
        to_fix = string_sel(fix_list,
                            names,
                            "Beam{0}".format(ibeam + 1))
        beam_param.set_fixed(to_fix)

      if beam_param.num_free() > 0:
        beam_params.append(beam_param)

    # Parameterise unique Crystals
    xl_ori_params = []
    xl_uc_params = []
    for icrystal, crystal in enumerate(experiments.crystals()):
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

      if options.scan_varying:
        if (not options.crystal.orientation.force_static or
            not options.crystal.unit_cell.force_static):
          # If a crystal is scan-varying, then it must always be found alongside
          # the same Scan and Goniometer in any Experiments in which it appears
          if [goniometer, scan].count(None) != 0:
            raise Sorry('A scan-varying crystal model cannot be created because '
                        'a scan or goniometer model is missing')
          if not all(g is goniometer and s is scan for (g, s) in assoc_models):
            raise Sorry('A single scan-varying crystal model cannot be refined '
                        'when associated with more than one scan or goniometer')
        sweep_range_deg = scan.get_oscillation_range(deg=True)
        array_range = scan.get_array_range()

        # orientation parameterisation
        if not options.crystal.orientation.force_static:
          n_intervals = options.crystal.orientation.smoother.absolute_num_intervals
          if n_intervals is None:
            deg_per_interval = options.crystal.orientation.smoother.interval_width_degrees
            if deg_per_interval is libtbx.Auto and analysis is not None:
              intervals = [analysis[i]['interval_width'] for i in exp_ids]
              deg_per_interval = min(intervals)
            elif deg_per_interval is None:
              deg_per_interval = 36.0
            n_intervals = max(int(
              abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)

          xl_ori_param = par.ScanVaryingCrystalOrientationParameterisation(
              crystal,
              array_range,
              n_intervals,
              experiment_ids=exp_ids)
        else: # force model to be static
          xl_ori_param = par.CrystalOrientationParameterisation(crystal,
                                                          experiment_ids=exp_ids)

        # unit cell parameterisation
        if not options.crystal.unit_cell.force_static:
          n_intervals = options.crystal.unit_cell.smoother.absolute_num_intervals
          if n_intervals is None:
            deg_per_interval = options.crystal.unit_cell.smoother.interval_width_degrees
            if deg_per_interval is libtbx.Auto and analysis is not None:
              intervals = [analysis[i]['interval_width'] for i in exp_ids]
              deg_per_interval = min(intervals)
            elif deg_per_interval is None:
              deg_per_interval = 36.0
            n_intervals = max(int(
              abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)

          set_errors = options.crystal.unit_cell.set_scan_varying_errors
          xl_uc_param = par.ScanVaryingCrystalUnitCellParameterisation(
              crystal,
              array_range,
              n_intervals,
              experiment_ids=exp_ids,
              set_state_uncertainties=set_errors)
        else: # force model to be static
          xl_uc_param = par.CrystalUnitCellParameterisation(crystal,
                                                          experiment_ids=exp_ids)
      else: # all models scan-static
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

      ori_fix_list = []
      if options.crystal.orientation.fix_list:
        ori_fix_list.extend(options.crystal.orientation.fix_list)

      cell_fix_list = []
      if options.crystal.unit_cell.fix_list:
        cell_fix_list.extend(options.crystal.unit_cell.fix_list)

      if options.crystal.fix:
        if options.crystal.fix == "all":
          xl_ori_param.set_fixed([True] * num_ori)
          xl_uc_param.set_fixed([True] * num_uc)
        elif options.crystal.fix == "cell":
          xl_uc_param.set_fixed([True] * num_uc)
        elif options.crystal.fix == "orientation":
          xl_ori_param.set_fixed([True] * num_ori)
        else: # can only get here if refinement.phil is broken
          raise RuntimeError("crystal.fix value not recognised")

      if cell_fix_list:
        names = filter_parameter_names(xl_uc_param)
        assert len(names) == num_uc
        to_fix = string_sel(cell_fix_list,
                            names,
                            "Crystal{0}".format(icrystal + 1))
        xl_uc_param.set_fixed(to_fix)

      if ori_fix_list:
        names = filter_parameter_names(xl_ori_param)
        assert len(names) == num_ori
        to_fix = string_sel(ori_fix_list,
                            names,
                            "Crystal{0}".format(icrystal + 1))
        xl_ori_param.set_fixed(to_fix)

      if xl_ori_param.num_free() > 0:
        xl_ori_params.append(xl_ori_param)
      if xl_uc_param.num_free() > 0:
        xl_uc_params.append(xl_uc_param)

    # Parameterise unique Detectors
    det_params = []
    for idetector, detector in enumerate(experiments.detectors()):
      # keep associated gonio and scan in case we are scan-varying
      exp_ids = experiments.indices(detector)
      assoc_models = [(experiments[i].goniometer, experiments[i].scan) \
                      for i in exp_ids]
      goniometer, scan = assoc_models[0]

      if options.detector.panels == "automatic":
        if len(detector) > 1:
          if options.scan_varying and not options.detector.force_static:
            raise Sorry('Scan-varying multiple panel detectors are not '
                        'currently supported')
          try:
            h = detector.hierarchy()
            det_param = par.DetectorParameterisationHierarchical(detector,
                experiment_ids=exp_ids, level=options.detector.hierarchy_level)
          except AttributeError:
            det_param = par.DetectorParameterisationMultiPanel(detector, beam,
                                                        experiment_ids=exp_ids)
        elif options.scan_varying and not options.detector.force_static:
          # If a detector is scan-varying, then it must always be found alongside
          # the same Scan and Goniometer in any Experiments in which it appears
          if [goniometer, scan].count(None) != 0:
            raise Sorry('A scan-varying detector model cannot be created '
                        'because a scan or goniometer model is missing')
          if not all(g is goniometer and s is scan for (g, s) in assoc_models):
            raise Sorry('A single scan-varying detector model cannot be '
              'refined when associated with more than one scan or goniometer')
          sweep_range_deg = scan.get_oscillation_range(deg=True)
          array_range = scan.get_array_range()

          n_intervals = options.detector.smoother.absolute_num_intervals
          if n_intervals is None:
            deg_per_interval = options.detector.smoother.interval_width_degrees
            if deg_per_interval is libtbx.Auto and analysis is not None:
              intervals = [analysis[i]['interval_width'] for i in exp_ids]
              deg_per_interval = min(intervals)
            elif deg_per_interval is None:
              deg_per_interval = 36.0
            n_intervals = max(int(
              abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)

          det_param = par.ScanVaryingDetectorParameterisationSinglePanel(
              detector,
              array_range,
              n_intervals,
              experiment_ids=exp_ids)
        else:
          det_param = par.DetectorParameterisationSinglePanel(detector,
                                                        experiment_ids=exp_ids)
      elif options.detector.panels == "single":
        if options.scan_varying and not options.detector.force_static:
          # If a detector is scan-varying, then it must always be found alongside
          # the same Scan and Goniometer in any Experiments in which it appears
          if [goniometer, scan].count(None) != 0:
            raise Sorry('A scan-varying detector model cannot be created '
                        'because a scan or goniometer model is missing')
          if not all(g is goniometer and s is scan for (g, s) in assoc_models):
            raise Sorry('A single scan-varying detector model cannot be '
              'refined when associated with more than one scan or goniometer')
          sweep_range_deg = scan.get_oscillation_range(deg=True)
          array_range = scan.get_array_range()

          n_intervals = options.detector.smoother.absolute_num_intervals
          if n_intervals is None:
            deg_per_interval = options.detector.smoother.interval_width_degrees
            if deg_per_interval is libtbx.Auto and analysis is not None:
              intervals = [analysis[i]['interval_width'] for i in exp_ids]
              deg_per_interval = min(intervals)
              for i in exp_ids:
                logger.debug(('Detector interval_width_degrees for experiment id'
                    ' {0} set to {1:.1f}').format(i, deg_per_interval))
            elif deg_per_interval is None:
              deg_per_interval = 36.0
            n_intervals = max(int(
              abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)

          det_param = par.ScanVaryingDetectorParameterisationSinglePanel(
              detector,
              array_range,
              n_intervals,
              experiment_ids=exp_ids)
        else:
          det_param = par.DetectorParameterisationSinglePanel(detector,
                                                        experiment_ids=exp_ids)
      elif options.detector.panels == "multiple":
        if options.scan_varying and not options.detector.force_static:
          raise Sorry('Scan-varying multiple panel detectors are not '
                      'currently supported')
        det_param = par.DetectorParameterisationMultiPanel(detector, beam,
                                                        experiment_ids=exp_ids)
      elif options.detector.panels == "hierarchical":
        if options.scan_varying and not options.detector.force_static:
          raise Sorry('Scan-varying hierarchical detectors are not '
                      'currently supported')
        det_param = par.DetectorParameterisationHierarchical(detector, beam,
                experiment_ids=exp_ids, level=options.detector.hierarchy_level)
      else: # can only get here if refinement.phil is broken
        raise RuntimeError("detector.panels value not recognised")

      # get number of fixable units, either parameters or parameter sets in
      # the scan-varying case
      try:
        num_det = det_param.num_sets()
      except AttributeError:
        num_det = det_param.num_total()

      fix_list = []
      if options.detector.fix_list:
        fix_list.extend(options.detector.fix_list)

      if options.detector.fix:
        if options.detector.fix == "all":
          det_param.set_fixed([True] * num_det)
        elif options.detector.fix == "position":
          fix_list.extend(['Dist', 'Shift1', 'Shift2'])
        elif options.detector.fix == "orientation":
          fix_list.extend(['Tau'])
        else: # can only get here if refinement.phil is broken
          raise RuntimeError("detector.fix value not recognised")

      if fix_list:
        names = filter_parameter_names(det_param)
        assert len(names) == num_det
        to_fix = string_sel(fix_list,
                            names,
                            "Detector{0}".format(idetector + 1))
        det_param.set_fixed(to_fix)

      if det_param.num_free() > 0:
        det_params.append(det_param)

    # Parameterise unique Goniometer setting matrices
    gon_params = []
    for igoniometer, goniometer in enumerate(experiments.goniometers()):
      if goniometer is None: continue
      # A Goniometer is parameterised with reference to the beam axis.
      # Use the first Beam this Goniometer is associated with.
      exp_ids = experiments.indices(goniometer)
      assoc_models = [(experiments[i].beam, experiments[i].scan) \
                      for i in exp_ids]
      beam, scan = assoc_models[0]

      if options.scan_varying:
        if not options.goniometer.force_static:
          # If a goniometer is scan-varying, then it must always be found
          # alongside the same Scan in any Experiments in which it appears
          if not scan:
            raise Sorry('A scan-varying goniometer model cannot be created '
                        'because a scan model is missing')
          if not all(s is scan for (g, s) in assoc_models):
            raise Sorry('A single scan-varying goniometer model cannot be '
                        'refined when associated with more than one scan')
        sweep_range_deg = scan.get_oscillation_range(deg=True)
        array_range = scan.get_array_range()

        if not options.goniometer.force_static:
          n_intervals = options.goniometer.smoother.absolute_num_intervals
          if n_intervals is None:
            deg_per_interval = options.goniometer.smoother.interval_width_degrees
            if deg_per_interval is libtbx.Auto and analysis is not None:
              intervals = [analysis[i]['interval_width'] for i in exp_ids]
              deg_per_interval = min(intervals)
            elif deg_per_interval is None:
              deg_per_interval = 36.0
            n_intervals = max(int(
              abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1)

          gon_param = par.ScanVaryingGoniometerParameterisation(
                                              goniometer,
                                              array_range,
                                              n_intervals,
                                              beam=beam,
                                              experiment_ids=exp_ids)
        else: # force model to be static
          gon_param = par.GoniometerParameterisation(goniometer, beam,
                                                       experiment_ids=exp_ids)
      else:
        # Parameterise scan static goniometer
        gon_param = par.GoniometerParameterisation(goniometer, beam,
                                                       experiment_ids=exp_ids)

      # get number of fixable units, either parameters or parameter sets in
      # the scan-varying case
      try:
        num_gon = gon_param.num_sets()
      except AttributeError:
        num_gon = gon_param.num_total()

      fix_list = []
      if options.goniometer.fix_list:
        fix_list.extend(options.goniometer.fix_list)

      if options.goniometer.fix:
        if "all" in options.goniometer.fix:
          gon_param.set_fixed([True] * num_gon)
        if "in_beam_plane" in options.goniometer.fix:
          fix_list.append('Gamma1')
        if "out_beam_plane" in options.goniometer.fix:
          fix_list.append('Gamma2')

      if fix_list:
        names = filter_parameter_names(gon_param)
        assert len(names) == num_gon
        to_fix = string_sel(fix_list,
                            names,
                            "Goniometer{0}".format(igoniometer + 1))
        gon_param.set_fixed(to_fix)

      if gon_param.num_free() > 0:
        gon_params.append(gon_param)

    # Parameter auto reduction options
    def model_nparam_minus_nref(p, reflections):
      cutoff = options.auto_reduction.min_nref_per_parameter * p.num_free()

      #Replaced Python code
      '''
      exp_ids = p.get_experiment_ids()
      # Do we have enough reflections to support this parameterisation?
      nparam = p.num_free()
      cutoff = options.auto_reduction.min_nref_per_parameter * nparam
      isel = flex.size_t()
      for exp_id in exp_ids:
        isel.extend((reflections['id'] == exp_id).iselection())
      nref = len(isel)
      return nref - cutoff
      '''
      return mnmn(reflections["id"],p.get_experiment_ids()).result - cutoff

    def unit_cell_nparam_minus_nref(p, reflections):
      '''Special version of model_nparam_minus_nref for crystal unit cell
      parameterisations. In some cases certain parameters of a unit cell
      parameterisation may affect only some subset of the total number of
      reflections. For example, for an orthorhombic cell the g_param_0 parameter
      has no effect on predictions in the plane (0,k,l). Here, take the number
      of affected reflections for each parameter into account.'''

      F_dbdp=flex.mat3_double( p.get_ds_dp() )
      min_nref = options.auto_reduction.min_nref_per_parameter
      # if no free parameters, do as model_nparam_minus_nref
      if len(F_dbdp) == 0:
        exp_ids = p.get_experiment_ids()
        isel = flex.size_t()
        for exp_id in exp_ids:
          isel.extend((reflections['id'] == exp_id).iselection())
        return len(isel)

      #Replaced Python code
      '''
      exp_ids = p.get_experiment_ids()
      isel = flex.size_t()
      for exp_id in exp_ids:
        isel.extend((reflections['id'] == exp_id).iselection())
      ref = reflections.select(isel)
      h = ref['miller_index'].as_vec3_double()
      dB_dp = p.get_ds_dp()
      # if no free parameters, do as model_nparam_minus_nref
      if len(dB_dp) == 0: return len(isel)
      nref_each_param = []
      min_nref = options.auto_reduction.min_nref_per_parameter
      for der in dB_dp:
        der_mat = flex.mat3_double(len(h), der.elems)
        tst = (der_mat * h).norms()
        nref_each_param.append((tst > 0.0).count(True))

      return min([nref - min_nref for nref in nref_each_param])
      '''
      return ucnmn(reflections["id"], reflections["miller_index"], p.get_experiment_ids(), F_dbdp).result - min_nref

    # In the scan-varying case we can't calculate dB_dp before composing the
    # model, so revert to the original function
    if options.scan_varying:
      unit_cell_nparam_minus_nref = model_nparam_minus_nref

    def panel_gp_nparam_minus_nref(p, pnl_ids, group, reflections, verbose=False):
      """
      :param p: ModelParameterisation; parameters in model
      :param pnl_ids: panel IDs
      :param group: group ID
      :panel reflections: flex table of reflections
      :panel verbose:
      :return: returns surplus {int}
      """
      exp_ids = p.get_experiment_ids() #Experiments parameterised by this ModelParameterisation
      # Do we have enough reflections to support this parameterisation?
      gp_params = [gp == group for gp in p.get_param_panel_groups()] #select the group ids for each param that matches the arg 'group'
      fixlist = p.get_fixed() # Get the fixed parameters; list says yes or no over all
      free_gp_params = [a and not b for a,b in zip(gp_params, fixlist)] #Free params is the total less the fixed
      nparam = free_gp_params.count(True)
      cutoff = options.auto_reduction.min_nref_per_parameter * nparam
      isel = flex.size_t()
      #Use Boost.Python extension module to replace below code
      surplus = pgnmn(reflections["id"], reflections["panel"], pnl_ids, exp_ids, cutoff).result

      #Replaced Python code
      '''
      for exp_id in exp_ids:
        sub_expID = (reflections['id'] == exp_id).iselection()
        sub_panels_expID = reflections['panel'].select(sub_expID)
        for pnl in pnl_ids:
          isel.extend(sub_expID.select(sub_panels_expID == pnl))
      nref = len(isel)
      surplus = nref - cutoff
      '''
      if surplus < 0 and verbose:
        logger.warning('{0} reflections on panels {1} with a cutoff of {2}'.format(nref, pnl_ids, cutoff))
      return surplus

    def weak_parameterisation_search(beam_params, xl_ori_params, xl_uc_params,
        det_params, gon_params, reflections):
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
          name = 'Beam{0}'.format(i + 1)
      for i, p in enumerate(xl_ori_params):
        net_nref = model_nparam_minus_nref(p, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          name = 'Crystal{0} orientation'.format(i + 1)
      for i, p in enumerate(xl_uc_params):
        net_nref = unit_cell_nparam_minus_nref(p, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          name = 'Crystal{0} unit cell'.format(i + 1)
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
              name = 'Detector{0}PanelGroup{1}'.format(i + 1, pnl_gp + 1)
        except Exception:
          net_nref = model_nparam_minus_nref(p, reflections)
          if net_nref < nref_deficit:
            nref_deficit = net_nref
            weak = p
            panels = None
            pnl_gp = None
            name = 'Detector{0}'.format(i + 1)
      for i, p in enumerate(gon_params):
        net_nref = model_nparam_minus_nref(p, reflections)
        if net_nref < nref_deficit:
          nref_deficit = net_nref
          weak = p
          name = 'Goniometer{0}'.format(i + 1)
      return {'parameterisation':weak,
              'panels':panels,
              'panel_group_id':pnl_gp,
              'name':name}

    # As a special case for detector metrology, try reducing the number of
    # detector parameters if there are too few for some panel group. If this is
    # unsuccessful, fail outright.
    if options.auto_reduction.detector_reduce:
      reduce_list = options.auto_reduction.detector_reduce_list
      for i, dp in enumerate(det_params):
        to_fix = flex.bool(dp.get_fixed())
        try: # test for hierarchical detector parameterisation
          pnl_groups = dp.get_panel_ids_by_group()
          for igp, gp in enumerate(pnl_groups):
            surplus = panel_gp_nparam_minus_nref(dp, gp, igp, reflections, verbose=True)
            if surplus < 0:
              msg = ('Require {0} more reflections to parameterise Detector{1} '
                     'panel group {2}').format(-1*surplus, i + 1, igp + 1)
              logger.warning(msg + '\nAttempting reduction of non-essential parameters')
              names = filter_parameter_names(dp)
              prefix = 'Group{0}'.format(igp + 1)
              reduce_this_group = [prefix + e for e in reduce_list]
              to_fix |= flex.bool(string_sel(reduce_this_group, names))
              # try again, and fail if still unsuccessful
              surplus = panel_gp_nparam_minus_nref(dp, gp, igp, reflections, verbose=True)
              if surplus < 0:
                msg = msg.format(-1*surplus, i + 1, igp + 1)
                raise Sorry(msg + '\nFailing.')
        except AttributeError:
          if model_nparam_minus_nref(dp, reflections) < 0:
            mdl = 'Detector{0}'.format(i + 1)
            msg = failmsg.format(mdl)
            raise Sorry(msg)
        dp.set_fixed(to_fix)

    if options.auto_reduction.action == 'fail':
      failmsg = 'Too few reflections to parameterise {0}'
      failmsg += '\nTry modifying refinement.parameterisation.auto_reduction options'
      for i, bp in enumerate(beam_params):
        if model_nparam_minus_nref(bp, reflections) < 0:
          mdl = 'Beam{0}'.format(i + 1)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

      for i, xlo in enumerate(xl_ori_params):
        if model_nparam_minus_nref(xlo, reflections) < 0:
          mdl = 'Crystal{0} orientation'.format(i + 1)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

      for i, xluc in enumerate(xl_uc_params):
        if unit_cell_nparam_minus_nref(xluc, reflections) < 0:
          mdl = 'Crystal{0} unit cell'.format(i + 1)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

      for i, dp in enumerate(det_params):
        try: # test for hierarchical detector parameterisation
          pnl_groups = dp.get_panel_ids_by_group()
          for igp, gp in enumerate(pnl_groups):
            if panel_gp_nparam_minus_nref(dp, gp, igp, reflections) < 0:
              msg = 'Too few reflections to parameterise Detector{0} panel group {1}'
              msg = msg.format(i + 1, igp + 1)
              msg += '\nTry modifying refinement.parameterisation.auto_reduction options'
              raise Sorry(msg)
        except AttributeError:
          if model_nparam_minus_nref(dp, reflections) < 0:
            mdl = 'Detector{0}'.format(i + 1)
            msg = failmsg.format(mdl)
            raise Sorry(msg)

      for i, gonp in enumerate(gon_params):
        if model_nparam_minus_nref(gonp, reflections) < 0:
          mdl = 'Goniometer{0}'.format(i + 1)
          msg = failmsg.format(mdl)
          raise Sorry(msg)

    elif options.auto_reduction.action == 'fix':
      warnmsg = 'Too few reflections to parameterise {0}'
      tmp = []
      for i, bp in enumerate(beam_params):
        if model_nparam_minus_nref(bp, reflections) >= 0:
          tmp.append(bp)
        else:
          mdl = 'Beam{0}'.format(i + 1)
          msg = warnmsg.format(mdl)
          logger.warning(msg)
      beam_params = tmp

      tmp = []
      for i, xlo in enumerate(xl_ori_params):
        if model_nparam_minus_nref(xlo, reflections) >= 0:
          tmp.append(xlo)
        else:
          mdl = 'Crystal{0} orientation'.format(i + 1)
          msg = warnmsg.format(mdl)
          logger.warning(msg)
      xl_ori_params = tmp

      tmp = []
      for i, xluc in enumerate(xl_uc_params):
        if unit_cell_nparam_minus_nref(xluc, reflections) >= 0:
          tmp.append(xluc)
        else:
          mdl = 'Crystal{0} unit cell'.format(i + 1)
          msg = warnmsg.format(mdl)
          logger.warning(msg)
      xl_uc_params = tmp

      tmp = []
      for i, dp in enumerate(det_params):
        fixlist = dp.get_fixed()
        try: # test for hierarchical detector parameterisation
          pnl_groups = dp.get_panel_ids_by_group()
          for igp, gp in enumerate(pnl_groups):
            if panel_gp_nparam_minus_nref(dp, gp, igp, reflections) < 0:
              msg = 'Too few reflections to parameterise Detector{0}PanelGroup{1}'
              msg = msg.format(i + 1, igp + 1)
              logger.warning(msg)
              gp_params = [gp == igp for gp in dp.get_param_panel_groups()]
              for j, val in enumerate(gp_params):
                if val: fixlist[j] = True
          dp.set_fixed(fixlist)
          if dp.num_free() > 0:
            tmp.append(dp)
          else:
            msg = 'No parameters remain free for Detector{0}'.format(i + 1)
            logger.warning(msg)
        except AttributeError:
          if model_nparam_minus_nref(dp, reflections) >= 0:
            tmp.append(dp)
          else:
            mdl = 'Detector{0}'.format(i + 1)
            msg = warnmsg.format(mdl)
            logger.warning(msg)
      det_params = tmp

      tmp = []
      for i, gonp in enumerate(gon_params):
        if model_nparam_minus_nref(gonp, reflections) >= 0:
          tmp.append(gonp)
        else:
          mdl = 'Goniometer{0}'.format(i + 1)
          msg = warnmsg.format(mdl)
          logger.warning(msg)
      gon_params = tmp

    elif options.auto_reduction.action == 'remove':
      # if there is only one experiment, it should be multi-panel for remove to make sense
      if len(experiments) == 1:
        if not det_params[-1].is_multi_state():
          raise Sorry("For single experiment, single panel refinement "
            "auto_reduction.action=remove cannot be used as it could only "
            "remove all reflections from refinement")
      warnmsg = 'Too few reflections to parameterise {0}'
      warnmsg += '\nAssociated reflections will be removed from the Reflection Manager'
      while True:
        dat = weak_parameterisation_search(beam_params, xl_ori_params,
            xl_uc_params, det_params, gon_params, reflections)
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
        logger.warning(msg)

      # Strip out parameterisations with zero free parameters
      beam_params = [p for p in beam_params if p.num_free() > 0]
      xl_ori_params = [p for p in xl_ori_params if p.num_free() > 0]
      xl_uc_params = [p for p in xl_uc_params if p.num_free() > 0]
      det_params = [p for p in det_params if p.num_free() > 0]
      gon_params = [p for p in gon_params if p.num_free() > 0]

    # Now we have the final list of model parameterisations, build a restraints
    # parameterisation (if requested). Only unit cell restraints are supported
    # at the moment.
    if any([options.crystal.unit_cell.restraints.tie_to_target,
            options.crystal.unit_cell.restraints.tie_to_group]):
      restraints_param = cls.config_restraints(params, det_params, beam_params,
        xl_ori_params, xl_uc_params, gon_params)
    else:
      restraints_param = None

    # Prediction equation parameterisation
    if do_stills: # doing stills
      if options.sparse:
        if options.spherical_relp_model:
          from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import SphericalRelpStillsPredictionParameterisationSparse as StillsPredictionParameterisation
        else:
          from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import StillsPredictionParameterisationSparse as StillsPredictionParameterisation
      else:
        if options.spherical_relp_model:
          from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
            import SphericalRelpStillsPredictionParameterisation as StillsPredictionParameterisation
        else:
          from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
              import StillsPredictionParameterisation
      pred_param = StillsPredictionParameterisation(
          experiments,
          det_params, beam_params, xl_ori_params, xl_uc_params)

    else: # doing scans
      if options.scan_varying:
        if options.sparse:
          from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
            import ScanVaryingPredictionParameterisationSparse as PredParam
        else:
          from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters \
            import ScanVaryingPredictionParameterisation as PredParam
        pred_param = PredParam(
              experiments,
              det_params, beam_params, xl_ori_params, xl_uc_params, gon_params)
      else:
        if options.sparse:
          from dials.algorithms.refinement.parameterisation.prediction_parameters \
            import XYPhiPredictionParameterisationSparse as PredParam
        else:
          from dials.algorithms.refinement.parameterisation.prediction_parameters \
            import XYPhiPredictionParameterisation as PredParam
        pred_param = PredParam(
            experiments,
            det_params, beam_params, xl_ori_params, xl_uc_params, gon_params)

    # Parameter reporting
    param_reporter = par.ParameterReporter(det_params, beam_params,
        xl_ori_params, xl_uc_params, gon_params)

    return pred_param, param_reporter, restraints_param

  @staticmethod
  def config_restraints(params, det_params, beam_params,
        xl_ori_params, xl_uc_params, gon_params):
    """Given a set of user parameters plus model parameterisations, create a
    restraints plus a parameterisation of these restraints

    Params:
        params The input parameters
        det_params A list of detector parameterisations
        beam_params A list of beam parameterisations
        xl_ori_params A list of crystal orientation parameterisations
        xl_uc_params A list of crystal unit cell parameterisations
        gon_params A list of goniometer parameterisations

    Returns:
        A restraints parameterisation
    """

    from dials.algorithms.refinement.restraints import RestraintsParameterisation
    rp = RestraintsParameterisation(detector_parameterisations = det_params,
               beam_parameterisations = beam_params,
               xl_orientation_parameterisations = xl_ori_params,
               xl_unit_cell_parameterisations = xl_uc_params,
               goniometer_parameterisations = gon_params)

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
      if tie.id is None:
        # get one experiment id for each parameterisation to apply to all
        tie.id = [e.get_experiment_ids()[0] for e in xl_uc_params]
      for exp_id in tie.id:
        rp.add_restraints_to_target_xl_unit_cell(exp_id, tie.values, tie.sigmas)

    for tie in cell_r.tie_to_group:
      if len(tie.sigmas) != 6:
        raise Sorry("6 sigmas must be provided as the tie_to_group.sigmas. "
                    "Note that individual sigmas of 0.0 will remove "
                    "the restraint for the corresponding cell parameter.")
      if tie.id is None:
        rp.add_restraints_to_group_xl_unit_cell(tie.target, "all", tie.sigmas)
      else:
        rp.add_restraints_to_group_xl_unit_cell(tie.target, tie.id, tie.sigmas)

    return rp

  @staticmethod
  def config_refinery(params, target, pred_param, constraints_manager,
    verbosity):
    """Given a set of parameters, a target class, a prediction
    parameterisation class and a constraints_manager (which could be None),
    build a refinery

    Params:
        params The input parameters

    Returns:
        The refinery instance
    """

    # Shorten parameter path
    options = params.refinement.refinery

    if options.engine == "SimpleLBFGS":
      from dials.algorithms.refinement.engine import SimpleLBFGS as refinery
    elif options.engine == "LBFGScurvs":
      from dials.algorithms.refinement.engine import LBFGScurvs as refinery
    elif options.engine == "GaussNewton":
      from dials.algorithms.refinement.engine import GaussNewtonIterations as refinery
    elif options.engine == "LevMar":
      from dials.algorithms.refinement.engine import LevenbergMarquardtIterations as refinery
    elif options.engine == "SparseLevMar":
      from dials.algorithms.refinement.sparse_engine import SparseLevenbergMarquardtIterations as refinery
    else:
      raise RuntimeError("Refinement engine " + options.engine +
                         " not recognised")

    logger.debug("Selected refinement engine type: %s", options.engine)

    engine = refinery(target = target,
            prediction_parameterisation = pred_param,
            constraints_manager=constraints_manager,
            log = options.log,
            verbosity = verbosity,
            tracking = options.journal,
            max_iterations = options.max_iterations)

    if params.refinement.mp.nproc > 1:
      nproc = params.refinement.mp.nproc
      try:
        engine.set_nproc(nproc)
      except NotImplementedError:
        logger.warning("Could not set nproc={0} for refinement engine of type {1}".format(
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

    # While a random subset of reflections is used, continue to
    # set random.seed to get consistent behaviour
    if options.random_seed is not None:
      import random
      random.seed(options.random_seed)
      flex.set_random_seed(options.random_seed)
      logger.debug("Random seed set to %d", options.random_seed)

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
      if options.weighting_strategy.override in ["stills", "external_deltapsi"]:
        msg = ('The "{0}" weighting strategy is not compatible with '
               'scan refinement').format(options.weighting_strategy.override)
        raise Sorry(msg)

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
        options.outlier.block_width=None
      else:
        colnames = ["x_resid", "y_resid", "phi_resid"]
      from dials.algorithms.refinement.outlier_detection import CentroidOutlierFactory
      outlier_detector = CentroidOutlierFactory.from_parameters_and_colnames(
        options, colnames, verbosity)

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
    elif options.weighting_strategy.override == "external_deltapsi":
      from dials.algorithms.refinement.weighting_strategies \
        import ExternalDelPsiWeightingStrategy
      weighting_strategy = ExternalDelPsiWeightingStrategy()
    elif options.weighting_strategy.override == "constant":
      from dials.algorithms.refinement.weighting_strategies \
        import ConstantWeightingStrategy
      weighting_strategy = ConstantWeightingStrategy(
        *options.weighting_strategy.constants, stills=do_stills)

    # calculate reflection block_width?
    if params.refinement.parameterisation.scan_varying:
      from dials.algorithms.refinement.reflection_manager import BlockCalculator
      block_calculator = BlockCalculator(experiments, reflections)
      if params.refinement.parameterisation.compose_model_per == "block":
        reflections = block_calculator.per_width(options.block_width, deg=True)
      elif params.refinement.parameterisation.compose_model_per == "image":
        reflections = block_calculator.per_image()

    return refman(reflections=reflections,
            experiments=experiments,
            nref_per_degree=options.reflections_per_degree,
            max_sample_size = options.maximum_sample_size,
            min_sample_size = options.minimum_sample_size,
            close_to_spindle_cutoff=options.close_to_spindle_cutoff,
            trim_scan_edges=options.trim_scan_edges,
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
    srm = params.refinement.parameterisation.spherical_relp_model

    if options.rmsd_cutoff == "fraction_of_bin_size":
      absolute_cutoffs = None
    elif options.rmsd_cutoff == "absolute":
      absolute_cutoffs = options.absolute_cutoffs
    else:
      raise RuntimeError("Target function rmsd_cutoff option" +
          options.rmsd_cutoff + " not recognised")

    # build managed reflection predictors
    from dials.algorithms.refinement.prediction import ExperimentsPredictor
    ref_predictor = ExperimentsPredictor(experiments, do_stills,
                                         spherical_relp=srm)

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

  def __init__(self, reflections, experiments,
               pred_param, param_reporter, refman, target, refinery,
               verbosity=0, copy_experiments=True):
    """
    Mandatory arguments:
      reflections - Input ReflectionList data
      experiments - a dxtbx ExperimentList object
      pred_param - An object derived from the PredictionParameterisation class
      param_reporter -A ParameterReporter object
      refman - A ReflectionManager object
      target - An object derived from the Target class
      refinery - An object derived from the Refinery class

    Optional arguments:
      verbosity - An integer verbosity level

    """

    # the experimental models
    self._experiments = experiments
    self.copy_experiments = copy_experiments

    # refinement module main objects
    self._pred_param = pred_param
    self._refman = refman
    self._target = target
    self._refinery = refinery

    # parameter reporter
    self._param_report = param_reporter

    self._verbosity = verbosity
    if verbosity == 0:
      logger.disabled = True

    return

  def get_experiments(self):

    if self.copy_experiments:
      from copy import deepcopy
      return deepcopy(self._experiments)
    else:
      return self._experiments

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

    corrmats = self._refinery.get_correlation_matrix_for_step(step)
    if corrmats is None: return None, None

    all_labels = self._pred_param.get_param_names()
    from dials.algorithms.refinement.refinement_helpers import string_sel
    if col_select is None:
      col_select = range(len(all_labels))
    sel = string_sel(col_select, all_labels)
    labels = [e for e, s in zip(all_labels, sel) if s]
    num_cols = num_rows = len(labels)
    if num_cols == 0: return None, None

    for k, corrmat in corrmats.items():

      assert corrmat.is_square_matrix()

      from scitbx.array_family import flex
      idx = flex.bool(sel).iselection()
      sub_corrmat = flex.double(flex.grid(num_cols, num_cols))

      for (i, x) in enumerate(idx):
        for (j, y) in enumerate(idx):
          sub_corrmat[i,j] = corrmat[x, y]

      corrmats[k] = sub_corrmat

    return (corrmats, labels)

  def print_step_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    logger.info("\nRefinement steps:")

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
    logger.info(st.format())
    logger.info(self._refinery.history.reason_for_termination)

    return

  def print_out_of_sample_rmsd_table(self):
    """print out-of-sample RSMDs per step, if these were tracked"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    # check if it makes sense to proceed
    if not "out_of_sample_rmsd" in self._refinery.history: return
    nref = len(self.get_free_reflections())
    if nref < 10: return # don't do anything if very few refs

    logger.info("\nRMSDs for out-of-sample (free) reflections:")

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
    logger.info(st.format())

    return

  def print_exp_rmsd_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    logger.info("\nRMSDs by experiment:")

    header = ["Exp\nid", "Nref"]
    for (name, units) in zip(self._target.rmsd_names, self._target.rmsd_units):
      if name == "RMSD_X" or name == "RMSD_Y" and units == "mm":
        header.append(name + "\n(px)")
      elif name == "RMSD_Phi" and units == "rad":
        # will convert radians to images for reporting of scans
        header.append("RMSD_Z" + "\n(images)")
      elif units == "rad":
        # will convert other angles in radians to degrees (e.g. for
        # RMSD_DeltaPsi and RMSD_2theta)
        header.append(name + "\n(deg)")
      else: # skip other/unknown RMSDs
        pass

    rows = []
    for iexp, exp in enumerate(self._experiments):
      detector = exp.detector
      px_sizes = [p.get_pixel_size() for p in detector]
      it = iter(px_sizes)
      px_size = next(it)
      if not all(tst == px_size for tst in it):
        logger.info("The detector in experiment %d does not have the same pixel " + \
             "sizes on each panel. Skipping...", iexp)
        continue
      px_per_mm = [1./e for e in px_size]

      scan = exp.scan
      try:
        temp = scan.get_oscillation(deg=False)
        images_per_rad  = 1./abs(scan.get_oscillation(deg=False)[1])
      except (AttributeError, ZeroDivisionError):
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
        elif units == "rad":
          rmsds.append(rmsd * rad2deg)
      rows.append([str(iexp), str(num)] + ["%.5g" % r for r in rmsds])

    if len(rows) > 0:
      truncated = False
      max_rows = 100
      if self._verbosity < 3 and len(rows) > max_rows:
        rows = rows[0:max_rows]
        truncated = True
      st = simple_table(rows, header)
      logger.info(st.format())
      if truncated:
        logger.info("Table truncated to show the first %d experiments only", max_rows)
        logger.info("Re-run with verbosity >= 3 to show all experiments")

    return

  def print_panel_rmsd_table(self):
    """print useful output about refinement steps in the form of a simple table"""

    from libtbx.table_utils import simple_table
    from math import pi
    rad2deg = 180/pi

    if len(self._experiments.scans()) > 1:
      logger.warning('Multiple scans present. Only the first scan will be used '
         'to determine the image width for reporting RMSDs')
    scan = self._experiments.scans()[0]
    try:
      temp = scan.get_oscillation(deg=False)
      images_per_rad  = 1./abs(scan.get_oscillation(deg=False)[1])
    except AttributeError:
      images_per_rad = None

    for idetector, detector in enumerate(self._experiments.detectors()):
      if len(detector) == 1: continue
      logger.info("\nDetector {0} RMSDs by panel:".format(idetector + 1))

      header = ["Panel\nid", "Nref"]
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
      for ipanel, panel in enumerate(detector):

        px_size = panel.get_pixel_size()
        px_per_mm = [1./e for e in px_size]
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
        logger.info(st.format())

    return

  def run(self):
    """Run refinement"""

    ####################################
    # Do refinement and return history #
    ####################################

    if self._verbosity > 1:
      logger.debug("\nExperimental models before refinement:")
      for i, beam in enumerate(self._experiments.beams()):
        logger.debug(ordinal_number(i) + ' ' + str(beam))
      for i, detector in enumerate(self._experiments.detectors()):
        logger.debug(ordinal_number(i) + ' ' + str(detector))
      for i, goniometer in enumerate(self._experiments.goniometers()):
        if goniometer is None: continue
        logger.debug(ordinal_number(i) + ' ' + str(goniometer))
      for i, scan in enumerate(self._experiments.scans()):
        if scan is None: continue
        logger.debug(ordinal_number(i) + ' ' + str(scan))
      for i, crystal in enumerate(self._experiments.crystals()):
        logger.debug(ordinal_number(i) + ' ' + str(crystal))

    self._refinery.run()

    # These involve calculation, so skip them when verbosity is zero, even
    # though the logger is disabled
    if self._verbosity > 0:
      self.print_step_table()
      self.print_out_of_sample_rmsd_table()
      self.print_exp_rmsd_table()

    det_npanels = [len(d) for d in self._experiments.detectors()]
    if any(n > 1 for n in det_npanels):
      self.print_panel_rmsd_table()

    # write scan-varying states back to their models
    #FIXME tidy up
    from dials.algorithms.refinement.parameterisation import \
      ScanVaryingPredictionParameterisation
    if isinstance(self._pred_param, ScanVaryingPredictionParameterisation):
      for iexp, exp in enumerate(self._experiments):
        ar_range = exp.scan.get_array_range()
        obs_image_numbers = range(ar_range[0], ar_range[1]+1)

        # write scan-varying s0 vectors back to beam models
        s0_list = self._pred_param.get_varying_s0(obs_image_numbers, iexp)
        if s0_list is not None:
          exp.beam.set_s0_at_scan_points(s0_list)

        # write scan-varying setting rotation matrices back to goniometer models
        S_list = self._pred_param.get_varying_setting_rotation(
            obs_image_numbers, iexp)
        if S_list is not None:
          exp.goniometer.set_setting_rotation_at_scan_points(S_list)

        # write scan-varying crystal setting matrices back to crystal models
        A_list = self._pred_param.get_varying_UB(obs_image_numbers, iexp)
        if A_list is not None:
          exp.crystal.set_A_at_scan_points(A_list)

        # get state covariance matrices the whole range of images. We select
        # the first element of this at each image because crystal scan-varying
        # parameterisations are not multi-state
        state_cov_list = [self._pred_param.calculate_model_state_uncertainties(
          obs_image_number=t, experiment_id=iexp) for t in range(ar_range[0],
                                                            ar_range[1]+1)]
        if 'U_cov' in state_cov_list[0]:
          u_cov_list = [e['U_cov'] for e in state_cov_list]
        else:
          u_cov_list = None

        if 'B_cov' in state_cov_list[0]:
          b_cov_list = [e['B_cov'] for e in state_cov_list]
        else:
          b_cov_list = None

        # return these to the model parameterisations to be set in the models
        self._pred_param.set_model_state_uncertainties(
          u_cov_list, b_cov_list, iexp)

    if self._verbosity > 1:
      logger.debug("\nExperimental models after refinement:")
      for i, beam in enumerate(self._experiments.beams()):
        logger.debug(ordinal_number(i) + ' ' + str(beam))
      for i, detector in enumerate(self._experiments.detectors()):
        logger.debug(ordinal_number(i) + ' ' + str(detector))
      for i, goniometer in enumerate(self._experiments.goniometers()):
        if goniometer is None: continue
        logger.debug(ordinal_number(i) + ' ' + str(goniometer))
      for i, scan in enumerate(self._experiments.scans()):
        if scan is None: continue
        logger.debug(ordinal_number(i) + ' ' + str(scan))
      for i, crystal in enumerate(self._experiments.crystals()):
        logger.debug(ordinal_number(i) + ' ' + str(crystal))

      # Report on the refined parameters
      logger.debug(str(self._param_report))

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
    refinement and additionally set the used_in_refinement flag. Do not
    compose the derivatives of states of the model as this is expensive and
    they are not needed outside of a refinement run"""

    reflections = self.predict_for_reflection_table(self._refman.get_indexed(),
      skip_derivatives=True)
    reflections.sort('iobs')
    mask = self.selection_used_for_refinement()
    reflections.set_flags(mask, reflections.flags.used_in_refinement)
    return reflections

  def predict_for_reflection_table(self, reflections, skip_derivatives=False):
    """perform prediction for all reflections in the supplied table"""

    # delegate to the target object, which has access to the predictor
    return self._target.predict_for_reflection_table(reflections, skip_derivatives)
