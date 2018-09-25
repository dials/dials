from __future__ import absolute_import, division
from dials.algorithms.refinement.parameterisation.beam_parameters import BeamParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.crystal_parameters import CrystalOrientationParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.crystal_parameters import CrystalUnitCellParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.detector_parameters import DetectorParameterisationSinglePanel # import dependency
from dials.algorithms.refinement.parameterisation.detector_parameters import DetectorParameterisationMultiPanel # import dependency
from dials.algorithms.refinement.parameterisation.detector_parameters import DetectorParameterisationHierarchical # import dependency
from dials.algorithms.refinement.parameterisation.goniometer_parameters import GoniometerParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.prediction_parameters import PredictionParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.prediction_parameters import XYPhiPredictionParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters import ScanVaryingCrystalOrientationParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters import ScanVaryingCrystalUnitCellParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.scan_varying_beam_parameters import ScanVaryingBeamParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.scan_varying_detector_parameters import ScanVaryingDetectorParameterisationSinglePanel # import dependency
from dials.algorithms.refinement.parameterisation.scan_varying_goniometer_parameters import ScanVaryingGoniometerParameterisation # import dependency
from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters import ScanVaryingPredictionParameterisation  # import dependency
from dials.algorithms.refinement.parameterisation.parameter_report import ParameterReporter # import dependency

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
