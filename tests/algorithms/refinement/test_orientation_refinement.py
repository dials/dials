"""
Test refinement of beam, detector and crystal orientation parameters
using generated reflection positions from ideal geometry.

Control of the experimental model and choice of minimiser is done via
PHIL, which means we can do, for example:

cctbx.python tst_orientation_refinement.py \
"random_seed=3; engine=LBFGScurvs"
"""


from __future__ import annotations

import sys
from math import pi
from os.path import join

from cctbx.sgtbx import space_group, space_group_symbols

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from dxtbx.format.FormatISISSXD import FormatISISSXD

# We will set up a mock scan and a mock experiment list
from dxtbx.model import CrystalFactory, ScanFactory
from dxtbx.model.experiment_list import Experiment, ExperimentList
from libtbx.phil import parse
from libtbx.test_utils import approx_equal
from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge
from scitbx import matrix
from scitbx.array_family import flex

from dials.algorithms.refinement.parameterisation.beam_parameters import (
    BeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import (
    DetectorParameterisationHierarchical,
    DetectorParameterisationSinglePanel,
)
from dials.algorithms.refinement.parameterisation.parameter_report import (
    ParameterReporter,
)

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    LauePredictionParameterisation,
    XYPhiPredictionParameterisation,
)
from dials.algorithms.refinement.prediction.managed_predictors import (
    LaueExperimentsPredictor,
    ScansExperimentsPredictor,
    ScansRayPredictor,
)
from dials.algorithms.refinement.refiner import Refiner, RefinerFactory
from dials.algorithms.refinement.reflection_manager import (
    LaueReflectionManager,
    ReflectionManager,
)

# Imports for the target function
from dials.algorithms.refinement.target import (
    LaueLeastSquaresResidualWithRmsdCutoff,
    LeastSquaresPositionalResidualWithRmsdCutoff,
)

# Reflection prediction
from dials.algorithms.spot_prediction import (
    IndexGenerator,
    LaueReflectionPredictor,
    ray_intersection,
)
from dials.command_line.refine import phil_scope

# Get modules to build models and minimiser using PHIL
from . import geometry_phil, minimiser_phil, setup_geometry, setup_minimiser


def test(args=[]):

    #############################
    # Setup experimental models #
    #############################

    master_phil = parse(f"{geometry_phil}\n{minimiser_phil}")

    models = setup_geometry.Extract(master_phil, cmdline_args=args)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 180 degree sequence
    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 1800),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(1800)),
        deg=True,
    )
    sequence_range = myscan.get_oscillation_range(deg=False)
    im_width = myscan.get_oscillation(deg=False)[1]
    assert sequence_range == (0.0, pi)
    assert approx_equal(im_width, 0.1 * pi / 180.0)

    # Build an experiment list
    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )

    ###########################
    # Parameterise the models #
    ###########################

    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)

    # Fix beam to the X-Z plane (imgCIF geometry), fix wavelength
    s0_param.set_fixed([True, False, True])

    # Fix crystal parameters
    # xluc_param.set_fixed([True, True, True, True, True, True])

    ########################################################################
    # Link model parameterisations together into a parameterisation of the #
    # prediction equation                                                  #
    ########################################################################

    pred_param = XYPhiPredictionParameterisation(
        experiments, [det_param], [s0_param], [xlo_param], [xluc_param]
    )

    ################################
    # Apply known parameter shifts #
    ################################

    # shift detector by 1.0 mm each translation and 2 mrad each rotation
    det_p_vals = det_param.get_param_vals()
    p_vals = [a + b for a, b in zip(det_p_vals, [1.0, 1.0, 1.0, 2.0, 2.0, 2.0])]
    det_param.set_param_vals(p_vals)

    # shift beam by 2 mrad in free axis
    s0_p_vals = s0_param.get_param_vals()
    p_vals = list(s0_p_vals)

    p_vals[0] += 2.0
    s0_param.set_param_vals(p_vals)

    # rotate crystal a bit (=2 mrad each rotation)
    xlo_p_vals = xlo_param.get_param_vals()
    p_vals = [a + b for a, b in zip(xlo_p_vals, [2.0, 2.0, 2.0])]
    xlo_param.set_param_vals(p_vals)

    # change unit cell a bit (=0.1 Angstrom length upsets, 0.1 degree of
    # gamma angle)
    xluc_p_vals = xluc_param.get_param_vals()
    cell_params = mycrystal.get_unit_cell().parameters()
    cell_params = [a + b for a, b in zip(cell_params, [0.1, 0.1, 0.1, 0.0, 0.0, 0.1])]
    new_uc = unit_cell(cell_params)
    newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
    S = symmetrize_reduce_enlarge(mycrystal.get_space_group())
    S.set_orientation(orientation=newB)
    X = tuple([e * 1.0e5 for e in S.forward_independent_parameters()])
    xluc_param.set_param_vals(X)

    #############################
    # Generate some reflections #
    #############################

    print("Reflections will be generated with the following geometry:")
    print(mybeam)
    print(mydetector)
    print(mycrystal)
    print("Target values of parameters are")
    msg = "Parameters: " + "%.5f " * len(pred_param)
    print(msg % tuple(pred_param.get_param_vals()))
    print()

    # All indices in a 2.0 Angstrom sphere
    resolution = 2.0
    index_generator = IndexGenerator(
        mycrystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()

    # Predict rays within the sequence range
    ray_predictor = ScansRayPredictor(experiments, sequence_range)
    obs_refs = ray_predictor(indices)

    print("Total number of reflections excited", len(obs_refs))

    # Take only those rays that intersect the detector
    intersects = ray_intersection(mydetector, obs_refs)
    obs_refs = obs_refs.select(intersects)

    # Make a reflection predictor and re-predict for all these reflections. The
    # result is the same, but we gain also the flags and xyzcal.px columns
    ref_predictor = ScansExperimentsPredictor(experiments)
    obs_refs["id"] = flex.int(len(obs_refs), 0)
    obs_refs = ref_predictor(obs_refs)

    # Set 'observed' centroids from the predicted ones
    obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.0
    px_size = mydetector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

    print("Total number of observations made", len(obs_refs))

    ###############################
    # Undo known parameter shifts #
    ###############################

    s0_param.set_param_vals(s0_p_vals)
    det_param.set_param_vals(det_p_vals)
    xlo_param.set_param_vals(xlo_p_vals)
    xluc_param.set_param_vals(xluc_p_vals)

    print("Initial values of parameters are")
    msg = "Parameters: " + "%.5f " * len(pred_param)
    print(msg % tuple(pred_param.get_param_vals()))
    print()

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(obs_refs, experiments)

    ##############################
    # Set up the target function #
    ##############################

    # The current 'achieved' criterion compares RMSD against 1/3 the pixel size and
    # 1/3 the image width in radians. For the simulated data, these are just made up

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
        experiments, ref_predictor, refman, pred_param, restraints_parameterisation=None
    )

    ################################
    # Set up the refinement engine #
    ################################

    refiner = setup_minimiser.Extract(
        master_phil, mytarget, pred_param, cmdline_args=args
    ).refiner

    print("Prior to refinement the experimental model is:")
    print(mybeam)
    print(mydetector)
    print(mycrystal)

    refiner.run()

    print()
    print("Refinement has completed with the following geometry:")
    print(mybeam)
    print(mydetector)
    print(mycrystal)


def test_laue_refinement(dials_data):
    fmt = FormatISISSXD(
        join(dials_data("isis_sxd_example_data", pathlib=True), "sxd_nacl_run.nxs")
    )
    beam = fmt.get_beam()
    detector = fmt.get_detector()
    goniometer = fmt.get_goniometer()
    scan = fmt.get_scan()
    crystal = CrystalFactory.from_dict(
        {
            "__id__": "crystal",
            "real_space_a": (
                0.5681647125795644,
                -2.9735716012061135,
                -2.707784412005687,
            ),
            "real_space_b": (
                -2.4994848902125884,
                -2.3900344014694066,
                2.091613643314567,
            ),
            "real_space_c": (
                -1.2771711635863638,
                3.676428861690809,
                -1.226011051463438,
            ),
            "space_group_hall_symbol": " P 1",
            "B_covariance": (
                2.618491627225783e-13,
                -2.4190170785778272e-30,
                2.7961382012436816e-30,
                1.4283218313839273e-13,
                8.110824693143866e-15,
                2.7961382012436816e-30,
                -1.922218398881239e-13,
                -1.1641948761717081e-14,
                2.2832201114561855e-14,
                -2.419017078577827e-30,
                1.3543505986455804e-44,
                -8.081590630292518e-46,
                -4.202632560757537e-29,
                -5.437640708903305e-29,
                -8.081590630292518e-46,
                3.330706229067803e-30,
                5.621471188408899e-29,
                -6.599119546892406e-30,
                2.7961382012436816e-30,
                -8.08159063029252e-46,
                9.550033948814972e-46,
                5.487666450546843e-30,
                2.7096475027184553e-30,
                9.550033948814972e-46,
                -3.935814660390771e-30,
                -3.889472044173952e-30,
                7.798194512461942e-30,
                1.428321831383927e-13,
                -4.2026325607575364e-29,
                5.487666450546843e-30,
                7.789867544667339e-13,
                1.4101250207277487e-13,
                5.487666450546843e-30,
                -2.0005409484272627e-13,
                -2.021584892435437e-13,
                4.481019714719027e-14,
                8.110824693143867e-15,
                -5.437640708903304e-29,
                2.7096475027184553e-30,
                1.4101250207277487e-13,
                2.5553690436147e-13,
                2.7096475027184553e-30,
                -1.1167612085554417e-14,
                -1.8848015530742402e-13,
                2.2125950964841596e-14,
                2.7961382012436816e-30,
                -8.08159063029252e-46,
                9.550033948814972e-46,
                5.487666450546843e-30,
                2.7096475027184553e-30,
                9.550033948814972e-46,
                -3.935814660390771e-30,
                -3.889472044173952e-30,
                7.798194512461942e-30,
                -1.922218398881239e-13,
                3.330706229067804e-30,
                -3.93581466039077e-30,
                -2.000540948427263e-13,
                -1.1167612085554417e-14,
                -3.93581466039077e-30,
                2.7092227778026175e-13,
                1.6029668235488112e-14,
                -3.2138365634328507e-14,
                -1.1641948761717081e-14,
                5.621471188408898e-29,
                -3.889472044173952e-30,
                -2.021584892435437e-13,
                -1.88480155307424e-13,
                -3.889472044173952e-30,
                1.6029668235488112e-14,
                2.7054780216756276e-13,
                -3.175994945548343e-14,
                2.2832201114561858e-14,
                -6.599119546892407e-30,
                7.79819451246194e-30,
                4.4810197147190265e-14,
                2.2125950964841592e-14,
                7.79819451246194e-30,
                -3.2138365634328507e-14,
                -3.175994945548343e-14,
                6.36770905528953e-14,
            ),
        }
    )

    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            crystal=crystal,
            imageset=None,
        )
    )

    det_param = DetectorParameterisationHierarchical(detector)
    xlo_param = CrystalOrientationParameterisation(crystal)
    xluc_param = CrystalUnitCellParameterisation(crystal)

    pred_param = LauePredictionParameterisation(
        experiments,
        detector_parameterisations=[det_param],
        beam_parameterisations=[],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param],
    )

    # shift detector by 0.2 mm each translation and 2 mrad each rotation
    det_p_vals = det_param.get_param_vals()
    p_vals = [a + b for a, b in zip(det_p_vals, [2.0, 2.0, 2.0, 2.0, 2.0, 2.0])]
    det_param.set_param_vals(p_vals)

    # rotate crystal a bit (=2 mrad each rotation)
    xlo_p_vals = xlo_param.get_param_vals()
    p_vals = [a + b for a, b in zip(xlo_p_vals, [2.0, 2.0, 2.0])]
    xlo_param.set_param_vals(p_vals)

    reflection_predictor = LaueReflectionPredictor(experiments[0], 1.0)
    obs_refs = reflection_predictor.all_reflections_for_asu(0.0)

    # Set 'observed' centroids from the predicted ones
    obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]
    obs_refs["s0"] = obs_refs["s0_cal"]
    obs_refs["wavelength"] = obs_refs["wavelength_cal"]
    obs_refs["id"] = flex.int(len(obs_refs), 0)

    # Invent some variances for the centroid positions of the simulated data
    px_size = detector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
    var_z = flex.double(len(obs_refs), 0.0)
    obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_z)

    # Undo known parameter shifts
    det_param.set_param_vals(det_p_vals)
    xlo_param.set_param_vals(xlo_p_vals)

    refman = LaueReflectionManager(obs_refs, experiments, outlier_detector=None)

    # Redefine the reflection predictor to use the type expected by the Target class
    ref_predictor = LaueExperimentsPredictor(experiments)

    target = LaueLeastSquaresResidualWithRmsdCutoff(
        experiments, ref_predictor, refman, pred_param, restraints_parameterisation=None
    )

    params = phil_scope.extract()
    param_reporter = ParameterReporter(
        pred_param.get_detector_parameterisations(),
        pred_param.get_beam_parameterisations(),
        pred_param.get_crystal_orientation_parameterisations(),
        pred_param.get_crystal_unit_cell_parameterisations(),
        pred_param.get_goniometer_parameterisations(),
    )
    refinery = RefinerFactory.config_refinery(params, target, pred_param, None)
    refiner = Refiner(experiments, pred_param, param_reporter, refman, target, refinery)

    print("Prior to refinement the experimental model is:")
    print(detector)
    print(crystal)

    refiner.run()

    print()
    print("Refinement has completed with the following geometry:")
    print(detector)
    print(crystal)


if __name__ == "__main__":
    test(sys.argv[1:])
