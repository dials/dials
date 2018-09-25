from __future__ import absolute_import, division
import logging
logger = logging.getLogger(__name__)

from dials_refinement_helpers_ext import pgnmn_iter as pgnmn
from dials_refinement_helpers_ext import ucnmn_iter as ucnmn
from dials_refinement_helpers_ext import mnmn_iter as mnmn

import libtbx # for libtbx.Auto
from scitbx.array_family import flex

# PHIL
from libtbx.phil import parse
from dials.algorithms.refinement.restraints.restraints_parameterisation \
  import uc_phil_str as uc_restraints_phil_str
from dials.algorithms.refinement.constraints import phil_str as constr_phil_str
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
  import phil_str as sv_phil_str

format_data = {'uc_restraints_phil':uc_restraints_phil_str,
               'constr_phil':constr_phil_str,
               'sv_phil_str':sv_phil_str}

phil_str = '''
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

    block_width = 1.0
      .help = "Width of a reflection 'block' (in degrees) determining how fine-"
              "grained the model used for scan-varying prediction during"
              "refinement is. Currently only has any effect if the crystal"
              "parameterisation is set to use compose_model_per=block"
      .type = float(value_min = 0.)
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
'''%format_data

phil_scope = parse(phil_str)


class ParameterisationFactory(object):

  @classmethod
  def from_parameters_and_experiments(cls, options, experiments,
      reflection_manager, do_stills=False):
    """Given a set of parameters, create a parameterisation from a set of
    experimental models.

    Params:
        options: The input parameters
        experiments: An ExperimentList object
        reflection_manager: A ReflectionManager object
        do_stills (bool)

    Returns:
        A tuple containing a prediction equation parameterisation and
        parameter reporter object.
    """

    # Shorten module paths
    import dials.algorithms.refinement.parameterisation as par

    # function to convert fix_lists into to_fix selections
    from dials.algorithms.refinement.refinement_helpers import string_sel

    # Get the working set of reflections
    reflections = reflection_manager.get_matches()

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
        ca = reflection_manager.get_centroid_analyser(debug=options.debug_centroid_analysis)
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
          obs = reflection_manager.get_obs()
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
          obs = reflection_manager.get_obs()
          isel=flex.size_t()
          for exp_id in exp_ids:
            isel.extend((obs['id'] == exp_id).iselection())
        # Now remove the selected reflections
        sel = flex.bool(len(obs), True)
        sel.set_selected(isel, False)
        reflection_manager.filter_obs(sel)
        reflections = reflection_manager.get_matches()
        logger.warning(msg)

      # Strip out parameterisations with zero free parameters
      beam_params = [p for p in beam_params if p.num_free() > 0]
      xl_ori_params = [p for p in xl_ori_params if p.num_free() > 0]
      xl_uc_params = [p for p in xl_uc_params if p.num_free() > 0]
      det_params = [p for p in det_params if p.num_free() > 0]
      gon_params = [p for p in gon_params if p.num_free() > 0]

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

    return pred_param, param_reporter
