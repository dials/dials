"""
Test refinement of beam, detector and crystal orientation parameters
using generated reflection positions from ideal geometry.

Control of the experimental model and choice of minimiser is done via
PHIL, which means we can do, for example:

cctbx.python tst_orientation_refinement.py \
"random_seed=3; engine=LBFGScurvs"
"""

# Python and cctbx imports
from __future__ import absolute_import, division, print_function
import sys


def test(args=[]):
    from math import pi
    from scitbx import matrix
    from scitbx.array_family import flex
    from libtbx.phil import parse
    from libtbx.test_utils import approx_equal

    # Get modules to build models and minimiser using PHIL
    import dials.test.algorithms.refinement.setup_geometry as setup_geometry
    import dials.test.algorithms.refinement.setup_minimiser as setup_minimiser

    # We will set up a mock scan and a mock experiment list
    from dxtbx.model import ScanFactory
    from dxtbx.model.experiment_list import ExperimentList, Experiment

    # Model parameterisations
    from dials.algorithms.refinement.parameterisation.detector_parameters import (
        DetectorParameterisationSinglePanel,
    )
    from dials.algorithms.refinement.parameterisation.beam_parameters import (
        BeamParameterisation,
    )
    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalOrientationParameterisation,
        CrystalUnitCellParameterisation,
    )

    # Symmetry constrained parameterisation for the unit cell
    from cctbx.uctbx import unit_cell
    from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge

    # Reflection prediction
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.refinement.prediction.managed_predictors import (
        ScansRayPredictor,
        ScansExperimentsPredictor,
    )
    from dials.algorithms.spot_prediction import ray_intersection
    from cctbx.sgtbx import space_group, space_group_symbols

    # Parameterisation of the prediction equation
    from dials.algorithms.refinement.parameterisation.prediction_parameters import (
        XYPhiPredictionParameterisation,
    )

    # Imports for the target function
    from dials.algorithms.refinement.target import (
        LeastSquaresPositionalResidualWithRmsdCutoff,
    )
    from dials.algorithms.refinement.reflection_manager import ReflectionManager

    #############################
    # Setup experimental models #
    #############################

    master_phil = parse(
        """
      include scope dials.test.algorithms.refinement.geometry_phil
      include scope dials.test.algorithms.refinement.minimiser_phil
      """,
        process_includes=True,
    )

    models = setup_geometry.Extract(master_phil, cmdline_args=args)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 180 degree sweep
    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 1800),
        exposure_times=0.1,
        oscillation=(0, 0.1),
        epochs=list(range(1800)),
        deg=True,
    )
    sweep_range = myscan.get_oscillation_range(deg=False)
    im_width = myscan.get_oscillation(deg=False)[1]
    assert sweep_range == (0.0, pi)
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

    # Predict rays within the sweep range
    ray_predictor = ScansRayPredictor(experiments, sweep_range)
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


if __name__ == "__main__":
    test(sys.argv[1:])
