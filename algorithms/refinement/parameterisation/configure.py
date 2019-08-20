from __future__ import absolute_import, division, print_function

import logging
import re

import libtbx  # for libtbx.Auto
from dials.algorithms.refinement import DialsRefineConfigError

# Function to convert fix_lists into to_fix selections
from dials.algorithms.refinement.refinement_helpers import string_sel

# Import model parameterisations
from .beam_parameters import BeamParameterisation
from .scan_varying_beam_parameters import ScanVaryingBeamParameterisation
from .crystal_parameters import CrystalOrientationParameterisation
from .scan_varying_crystal_parameters import (
    ScanVaryingCrystalOrientationParameterisation,
)
from .crystal_parameters import CrystalUnitCellParameterisation
from .scan_varying_crystal_parameters import ScanVaryingCrystalUnitCellParameterisation
from .detector_parameters import DetectorParameterisationHierarchical
from .detector_parameters import DetectorParameterisationMultiPanel
from .detector_parameters import DetectorParameterisationSinglePanel
from .scan_varying_detector_parameters import (
    ScanVaryingDetectorParameterisationSinglePanel,
)
from .goniometer_parameters import GoniometerParameterisation
from .scan_varying_goniometer_parameters import ScanVaryingGoniometerParameterisation

# Import parameterisations of the prediction equation
from .prediction_parameters import XYPhiPredictionParameterisation
from .prediction_parameters import XYPhiPredictionParameterisationSparse
from .scan_varying_prediction_parameters import ScanVaryingPredictionParameterisation
from .scan_varying_prediction_parameters import (
    ScanVaryingPredictionParameterisationSparse,
)
from .prediction_parameters_stills import StillsPredictionParameterisation
from .prediction_parameters_stills import StillsPredictionParameterisationSparse
from .prediction_parameters_stills import SphericalRelpStillsPredictionParameterisation
from .prediction_parameters_stills import (
    SphericalRelpStillsPredictionParameterisationSparse,
)

from .autoreduce import AutoReduce

# PHIL
from libtbx.phil import parse
from dials.algorithms.refinement.parameterisation.autoreduce import (
    phil_str as autoreduce_phil_str,
)
from dials.algorithms.refinement.restraints.restraints_parameterisation import (
    uc_phil_str as uc_restraints_phil_str,
)
from dials.algorithms.refinement.constraints import phil_str as constr_phil_str
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters import (
    phil_str as sv_phil_str,
)

format_data = {
    "autoreduce_phil": autoreduce_phil_str,
    "uc_restraints_phil": uc_restraints_phil_str,
    "constr_phil": constr_phil_str,
    "sv_phil": sv_phil_str,
}

logger = logging.getLogger(__name__)

phil_str = (
    """
    auto_reduction
      .help = "determine behaviour when there are too few reflections to"
              "reasonably produce a full parameterisation of the experiment list"
      .expert_level = 1
    {
      %(autoreduce_phil)s
    }

    scan_varying = False
      .help = "Allow models that are not forced to be static to vary during"
              "the scan, Auto will run one macrocycle with static then"
              "scan varying refinement for the crystal"
      .type = bool
      .short_caption = "Scan-varying refinement"

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

    set_scan_varying_errors = False
      .type = bool
      .help = "If scan-varying refinement is done, and if the estimated"
              "covariance of the model states have been calculated by the"
              "minimiser, choose whether to return this to the models or"
              "not. The default is not to, in order to keep the file size"
              "of the serialized model small. At the moment, this only has"
              "an effect for crystal unit cell (B matrix) errors."

    debug_centroid_analysis = False
      .help = "Set True to write out a file containing the reflections used"
              "for centroid analysis for automatic setting of the  scan-varying"
              "interval width. This can then be analysed with"
              "dev.dials.plot_centroid_analysis (requires dials_scratch repository)."
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
        .short_caption = "Fix beam parameters"

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

      %(sv_phil)s
    }

    crystal
      .help = "crystal parameters"
    {
      fix = all cell orientation
        .help = "Fix crystal parameters"
        .type = choice
        .short_caption = "Fix crystal parameters"

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

        %(sv_phil)s
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

        %(sv_phil)s
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
        .short_caption = "Fix detector parameters"

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

      %(sv_phil)s
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
        .short_caption = "Fix goniometer parameters"

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

      %(sv_phil)s
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
"""
    % format_data
)
phil_scope = parse(phil_str)


# A helper function for parameter fixing
def _filter_parameter_names(parameterisation):
    # scan-varying suffixes like '_sample1' should be removed from
    # the list of parameter names so that it is num_free in length
    pattern = re.compile(r"_sample[0-9]+$")
    names = [
        pattern.sub("", e) for e in parameterisation.get_param_names(only_free=False)
    ]
    filtered_names = []
    for name in names:
        if name not in filtered_names:
            filtered_names.append(name)
    return filtered_names


# Helper function to perform centroid analysis
def _centroid_analysis(options, experiments, reflection_manager):

    analysis = None
    if not options.scan_varying:
        return analysis

    tst = [
        options.beam.smoother,
        options.crystal.orientation.smoother,
        options.crystal.unit_cell.smoother,
        options.detector.smoother,
        options.goniometer.smoother,
    ]
    tst = [
        (e.absolute_num_intervals is None and e.interval_width_degrees is libtbx.Auto)
        for e in tst
    ]
    if any(tst):
        logger.info(
            "Doing centroid analysis to "
            "automatically determine scan-varying interval widths"
        )
        ca = reflection_manager.get_centroid_analyser(
            debug=options.debug_centroid_analysis
        )
        analysis = ca()
    if analysis is None:
        return analysis

    # Use the results of centroid analysis to suggest suitable interval widths
    # for each experiment. This will be the smallest of the proposed intervals
    # for each of the residuals in x, y and phi, as long as this is not smaller
    # than either the outlier rejection block width, or 9.0 degrees.
    for i, a in enumerate(analysis):
        intervals = (a.get("x_interval"), a.get("y_interval"), a.get("phi_interval"))
        try:
            min_interval = min(_f for _f in intervals if _f is not None)
        except ValueError:
            # empty list - analysis was unable to suggest a suitable interval
            # width. Default to the safest case
            phi_min, phi_max = experiments[i].scan.get_oscillation_range(deg=True)
            a["interval_width"] = abs(phi_max - phi_min)
            logger.info(
                "Exp id {0} suggested interval width could not be "
                "determined and will be reset to the scan width of "
                "{1:.1f} degrees".format(i, a["interval_width"])
            )
            continue
        min_interval = max(min_interval, 9.0)
        block_size = a.get("block_size")
        if block_size is not None:
            min_interval = max(min_interval, block_size)
        a["interval_width"] = min_interval
        logger.info(
            "Exp id {0} suggested interval width = {1:.1f} degrees".format(
                i, min_interval
            )
        )

    return analysis


# Helper function to choose the number of intervals for the smoother
def _set_n_intervals(smoother_params, analysis, scan, exp_ids):
    n_intervals = smoother_params.absolute_num_intervals
    if n_intervals is not None:
        return n_intervals

    deg_per_interval = smoother_params.interval_width_degrees
    if deg_per_interval is libtbx.Auto and analysis is not None:
        intervals = [analysis[i]["interval_width"] for i in exp_ids]
        deg_per_interval = min(intervals)
    if deg_per_interval is None:
        deg_per_interval = 36.0

    sweep_range_deg = scan.get_oscillation_range(deg=True)
    n_intervals = max(
        int(abs(sweep_range_deg[1] - sweep_range_deg[0]) / deg_per_interval), 1
    )
    return n_intervals


def _parameterise_beams(options, experiments, analysis):
    beam_params = []
    sv_beam = options.scan_varying and not options.beam.force_static
    for ibeam, beam in enumerate(experiments.beams()):
        # The Beam is parameterised with reference to a goniometer axis (or None).
        # Use the first (if any) Goniometers this Beam is associated with.
        exp_ids = experiments.indices(beam)
        assoc_models = [
            (experiments[i].goniometer, experiments[i].scan) for i in exp_ids
        ]
        goniometer, scan = assoc_models[0]

        if sv_beam:
            if not all((goniometer, scan)):
                raise DialsRefineConfigError(
                    "A scan-varying beam model cannot be created because "
                    "a scan or goniometer model is missing"
                )
            # If a beam is scan-varying, then it must always be found alongside
            # the same Scan and Goniometer in any Experiments in which it appears
            if not all(g is goniometer and s is scan for (g, s) in assoc_models):
                raise DialsRefineConfigError(
                    "A single scan-varying beam model cannot be refined "
                    "when associated with more than one scan or goniometer"
                )
            array_range = scan.get_array_range()
            n_intervals = _set_n_intervals(
                options.beam.smoother, analysis, scan, exp_ids
            )
            beam_param = ScanVaryingBeamParameterisation(
                beam,
                array_range,
                n_intervals,
                goniometer=goniometer,
                experiment_ids=exp_ids,
            )
        else:
            # Parameterise scan static beam, passing the goniometer
            beam_param = BeamParameterisation(beam, goniometer, experiment_ids=exp_ids)

        # Set the model identifier to name the parameterisation
        beam_param.model_identifier = "Beam{}".format(ibeam + 1)

        # get number of fixable units, either parameters or parameter sets in
        # the scan-varying case
        num_beam = getattr(beam_param, "num_sets", getattr(beam_param, "num_total"))()

        fix_list = []
        if options.beam.fix_list:
            fix_list.extend(options.beam.fix_list)

        if options.beam.fix:
            if "all" in options.beam.fix:
                beam_param.set_fixed([True] * num_beam)
            if "in_spindle_plane" in options.beam.fix:
                fix_list.append("Mu1")
            if "out_spindle_plane" in options.beam.fix:
                fix_list.append("Mu2")
            if "wavelength" in options.beam.fix:
                fix_list.append("nu")

        if fix_list:
            names = _filter_parameter_names(beam_param)
            assert len(names) == num_beam
            to_fix = string_sel(fix_list, names, beam_param.model_identifier)
            beam_param.set_fixed(to_fix)

        if beam_param.num_free() > 0:
            beam_params.append(beam_param)

    return beam_params


def _parameterise_crystals(options, experiments, analysis):
    xl_ori_params = []
    xl_uc_params = []
    sv_xl_ori = options.scan_varying and not options.crystal.orientation.force_static
    sv_xl_uc = options.scan_varying and not options.crystal.unit_cell.force_static
    for icrystal, crystal in enumerate(experiments.crystals()):
        # This crystal can only ever appear either in scans or in stills
        # (otherwise it requires a different crystal model)
        exp_ids = experiments.indices(crystal)
        assoc_models = [
            (experiments[i].goniometer, experiments[i].scan) for i in exp_ids
        ]
        goniometer, scan = assoc_models[0]
        if goniometer is None:
            # There should be no associated goniometer and scan models
            if any(g or s for (g, s) in assoc_models):
                raise DialsRefineConfigError(
                    "A crystal model appears in a mixture of scan and still "
                    "experiments, which is not supported"
                )

        if sv_xl_ori or sv_xl_uc:
            if not all((goniometer, scan)):
                raise DialsRefineConfigError(
                    "A scan-varying crystal model cannot be created because "
                    "a scan or goniometer model is missing"
                )
            # If a crystal is scan-varying, then it must always be found alongside
            # the same Scan and Goniometer in any Experiments in which it appears
            if not all(g is goniometer and s is scan for (g, s) in assoc_models):
                raise DialsRefineConfigError(
                    "A single scan-varying crystal model cannot be refined "
                    "when associated with more than one scan or goniometer"
                )
            array_range = scan.get_array_range()

        # orientation parameterisation
        if sv_xl_ori:
            n_intervals = _set_n_intervals(
                options.crystal.orientation.smoother, analysis, scan, exp_ids
            )
            xl_ori_param = ScanVaryingCrystalOrientationParameterisation(
                crystal, array_range, n_intervals, experiment_ids=exp_ids
            )
        else:  # force model to be static
            xl_ori_param = CrystalOrientationParameterisation(
                crystal, experiment_ids=exp_ids
            )

        # Set the model identifier to name the parameterisation
        xl_ori_param.model_identifier = "Crystal{}".format(icrystal + 1)

        # unit cell parameterisation
        if sv_xl_uc:
            n_intervals = _set_n_intervals(
                options.crystal.unit_cell.smoother, analysis, scan, exp_ids
            )
            set_errors = options.set_scan_varying_errors
            xl_uc_param = ScanVaryingCrystalUnitCellParameterisation(
                crystal,
                array_range,
                n_intervals,
                experiment_ids=exp_ids,
                set_state_uncertainties=set_errors,
            )
        else:  # force model to be static
            xl_uc_param = CrystalUnitCellParameterisation(
                crystal, experiment_ids=exp_ids
            )

        # Set the model identifier to name the parameterisation
        xl_uc_param.model_identifier = "Crystal{}".format(icrystal + 1)

        # get number of fixable units, either parameters or parameter sets in
        # the scan-varying case
        num_ori = getattr(
            xl_ori_param, "num_sets", getattr(xl_ori_param, "num_total")
        )()
        num_uc = getattr(xl_uc_param, "num_sets", getattr(xl_uc_param, "num_total"))()

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
            else:  # can only get here if refinement.phil is broken
                raise RuntimeError("crystal.fix value not recognised")

        if cell_fix_list:
            names = _filter_parameter_names(xl_uc_param)
            assert len(names) == num_uc
            to_fix = string_sel(cell_fix_list, names, xl_uc_param.model_identifier)
            xl_uc_param.set_fixed(to_fix)

        if ori_fix_list:
            names = _filter_parameter_names(xl_ori_param)
            assert len(names) == num_ori
            to_fix = string_sel(ori_fix_list, names, xl_ori_param.model_identifier)
            xl_ori_param.set_fixed(to_fix)

        if xl_ori_param.num_free() > 0:
            xl_ori_params.append(xl_ori_param)
        if xl_uc_param.num_free() > 0:
            xl_uc_params.append(xl_uc_param)

    return xl_ori_params, xl_uc_params


def _parameterise_detectors(options, experiments, analysis):
    det_params = []
    sv_det = options.scan_varying and not options.detector.force_static
    for idetector, detector in enumerate(experiments.detectors()):
        # keep associated gonio and scan in case we are scan-varying
        exp_ids = experiments.indices(detector)
        assoc_models = [
            (experiments[i].goniometer, experiments[i].scan) for i in exp_ids
        ]
        goniometer, scan = assoc_models[0]

        if sv_det:
            if not all((goniometer, scan)):
                raise DialsRefineConfigError(
                    "A scan-varying detector model cannot be created "
                    "because a scan or goniometer model is missing"
                )
            # If a detector is scan-varying, then it must always be found alongside
            # the same Scan and Goniometer in any Experiments in which it appears
            if not all(g is goniometer and s is scan for (g, s) in assoc_models):
                raise DialsRefineConfigError(
                    "A single scan-varying detector model cannot be "
                    "refined when associated with more than one scan or goniometer"
                )

            # Additional checks on whether a scan-varying parameterisation is allowed
            if options.detector.panels == "automatic" and len(detector) > 1:
                raise DialsRefineConfigError(
                    "Scan-varying multiple panel detectors are not "
                    "currently supported"
                )
            if options.detector.panels == "multiple":
                raise DialsRefineConfigError(
                    "Scan-varying multiple panel detectors are not "
                    "currently supported"
                )
            if options.detector.panels == "hierarchical":
                raise DialsRefineConfigError(
                    "Scan-varying hierarchical detectors are not " "currently supported"
                )

            array_range = scan.get_array_range()
            n_intervals = _set_n_intervals(
                options.detector.smoother, analysis, scan, exp_ids
            )
            det_param = ScanVaryingDetectorParameterisationSinglePanel(
                detector, array_range, n_intervals, experiment_ids=exp_ids
            )
        else:
            # Convert automatic into correct specific option
            if options.detector.panels == "automatic":
                if len(detector) > 1:
                    if hasattr(detector, "hierarchy"):
                        options.detector.panels = "hierarchical"
                    else:
                        options.detector.panels = "multiple"
                else:
                    options.detector.panels = "single"

            # Construct parameterisation based on panels choice
            if options.detector.panels == "single":
                if len(detector) > 1:
                    raise DialsRefineConfigError(
                        "A single panel parameterisation cannot be created "
                        "for a multiple panel detector"
                    )
                det_param = DetectorParameterisationSinglePanel(
                    detector, experiment_ids=exp_ids
                )
            elif options.detector.panels == "multiple":
                # Take first associated beam model
                beam = experiments[exp_ids[0]].beam
                det_param = DetectorParameterisationMultiPanel(
                    detector, beam, experiment_ids=exp_ids
                )
            elif options.detector.panels == "hierarchical":
                try:  # Use hierarchy in parameterisation if the detector has one
                    detector.hierarchy()
                    det_param = DetectorParameterisationHierarchical(
                        detector,
                        experiment_ids=exp_ids,
                        level=options.detector.hierarchy_level,
                    )
                except AttributeError:
                    raise DialsRefineConfigError(
                        "A hierarchical detector parameterisation cannot be "
                        "created for a detector without a hierarchy"
                    )

        # Set the model identifier to name the parameterisation
        det_param.model_identifier = "Detector{}".format(idetector + 1)

        # get number of fixable units, either parameters or parameter sets in
        # the scan-varying case
        num_det = getattr(det_param, "num_sets", getattr(det_param, "num_total"))()

        fix_list = []
        if options.detector.fix_list:
            fix_list.extend(options.detector.fix_list)

        if options.detector.fix:
            if options.detector.fix == "all":
                det_param.set_fixed([True] * num_det)
            elif options.detector.fix == "position":
                fix_list.extend(["Dist", "Shift1", "Shift2"])
            elif options.detector.fix == "orientation":
                fix_list.extend(["Tau"])
            else:  # can only get here if refinement.phil is broken
                raise RuntimeError("detector.fix value not recognised")

        if fix_list:
            names = _filter_parameter_names(det_param)
            assert len(names) == num_det
            to_fix = string_sel(fix_list, names, det_param.model_identifier)
            det_param.set_fixed(to_fix)

        if det_param.num_free() > 0:
            det_params.append(det_param)

    return det_params


def _parameterise_goniometers(options, experiments, analysis):
    gon_params = []
    sv_gon = options.scan_varying and not options.goniometer.force_static
    for igoniometer, goniometer in enumerate(experiments.goniometers()):
        if goniometer is None:
            continue
        # A Goniometer is parameterised with reference to the beam axis.
        # Use the first Beam this Goniometer is associated with.
        exp_ids = experiments.indices(goniometer)
        assoc_models = [(experiments[i].beam, experiments[i].scan) for i in exp_ids]
        beam, scan = assoc_models[0]

        if sv_gon:
            # If a goniometer is scan-varying, then it must always be found
            # alongside the same Scan in any Experiments in which it appears
            if not scan:
                raise DialsRefineConfigError(
                    "A scan-varying goniometer model cannot be created "
                    "because a scan model is missing"
                )
            if not all(s is scan for (g, s) in assoc_models):
                raise DialsRefineConfigError(
                    "A single scan-varying goniometer model cannot be "
                    "refined when associated with more than one scan"
                )
            array_range = scan.get_array_range()
            n_intervals = _set_n_intervals(
                options.goniometer.smoother, analysis, scan, exp_ids
            )
            gon_param = ScanVaryingGoniometerParameterisation(
                goniometer, array_range, n_intervals, beam=beam, experiment_ids=exp_ids
            )
        else:  # force model to be static
            gon_param = GoniometerParameterisation(
                goniometer, beam, experiment_ids=exp_ids
            )

        # Set the model identifier to name the parameterisation
        gon_param.model_identifier = "Goniometer{}".format(igoniometer + 1)

        # get number of fixable units, either parameters or parameter sets in
        # the scan-varying case
        num_gon = getattr(gon_param, "num_sets", getattr(gon_param, "num_total"))()

        fix_list = []
        if options.goniometer.fix_list:
            fix_list.extend(options.goniometer.fix_list)

        if options.goniometer.fix:
            if "all" in options.goniometer.fix:
                gon_param.set_fixed([True] * num_gon)
            if "in_beam_plane" in options.goniometer.fix:
                fix_list.append("Gamma1")
            if "out_beam_plane" in options.goniometer.fix:
                fix_list.append("Gamma2")

        if fix_list:
            names = _filter_parameter_names(gon_param)
            assert len(names) == num_gon
            to_fix = string_sel(fix_list, names, gon_param.model_identifier)
            gon_param.set_fixed(to_fix)

        if gon_param.num_free() > 0:
            gon_params.append(gon_param)

    return gon_params


def build_prediction_parameterisation(
    options, experiments, reflection_manager, do_stills=False
):
    """Given a set of parameters, create a parameterisation from a set of
    experimental models.

    Params:
        options: The input parameters
        experiments: An ExperimentList object
        reflection_manager: A ReflectionManager object
        do_stills (bool)

    Returns:
        A prediction equation parameterisation object
    """

    # If required, do full centroid analysis on the reflections (assumes
    # outlier-rejection has been done already) to determine suitable interval
    # widths for scan-varying refinement
    analysis = _centroid_analysis(options, experiments, reflection_manager)

    # Parameterise each unique model
    beam_params = _parameterise_beams(options, experiments, analysis)
    xl_ori_params, xl_uc_params = _parameterise_crystals(options, experiments, analysis)
    det_params = _parameterise_detectors(options, experiments, analysis)
    gon_params = _parameterise_goniometers(options, experiments, analysis)

    # Check for too many parameters and reduce if requested
    autoreduce = AutoReduce(
        options.auto_reduction,
        det_params,
        beam_params,
        xl_ori_params,
        xl_uc_params,
        gon_params,
        reflection_manager,
        scan_varying=options.scan_varying,
    )
    autoreduce()
    det_params = autoreduce.det_params
    beam_params = autoreduce.beam_params
    xl_ori_params = autoreduce.xl_ori_params
    xl_uc_params = autoreduce.xl_uc_params
    gon_params = autoreduce.gon_params

    # Build the prediction equation parameterisation
    if do_stills:  # doing stills
        if options.sparse:
            if options.spherical_relp_model:
                PredParam = SphericalRelpStillsPredictionParameterisationSparse
            else:
                PredParam = StillsPredictionParameterisationSparse
        else:
            if options.spherical_relp_model:
                PredParam = SphericalRelpStillsPredictionParameterisation
            else:
                PredParam = StillsPredictionParameterisation
        pred_param = PredParam(
            experiments, det_params, beam_params, xl_ori_params, xl_uc_params
        )

    else:  # doing scans
        if options.scan_varying:
            if options.sparse:
                PredParam = ScanVaryingPredictionParameterisationSparse
            else:
                PredParam = ScanVaryingPredictionParameterisation
            pred_param = PredParam(
                experiments,
                det_params,
                beam_params,
                xl_ori_params,
                xl_uc_params,
                gon_params,
            )
            pred_param.set_scan_varying_errors = options.set_scan_varying_errors

        else:
            if options.sparse:
                PredParam = XYPhiPredictionParameterisationSparse
            else:
                PredParam = XYPhiPredictionParameterisation
            pred_param = PredParam(
                experiments,
                det_params,
                beam_params,
                xl_ori_params,
                xl_uc_params,
                gon_params,
            )

    return pred_param
