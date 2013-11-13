#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Investigate the convergence radius of refinement by multiple refinement runs
doing random sampling of parameter space."""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get module to build experimental models
from setup_geometry import Extract

# We will set up a mock scan
from dxtbx.model.scan import scan_factory

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation
from dials.algorithms.refinement.refinement_helpers import random_param_shift

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints.parameter_reduction import \
    symmetrize_reduce_enlarge

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ReflectionPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    DetectorSpacePredictionParameterisation

# Imports for the target function
from dials.algorithms.refinement.target import \
    LeastSquaresPositionalResidualWithRmsdCutoff, ReflectionManager

# Import the refinement engine
from dials.algorithms.refinement.engine import SimpleLBFGS
from dials.algorithms.refinement.engine import LBFGScurvs
from dials.algorithms.refinement.engine import GaussNewtonIterations

# Import helper functions
from dials.algorithms.refinement.refinement_helpers import print_model_geometry

def setup_models(seed):

    #############################
    # Setup experimental models #
    #############################

    override = "geometry.parameters.random_seed=" + str(seed)
    print override

    master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    include scope dials.test.algorithms.refinement.minimiser_phil
    """, process_includes=True)

    models = Extract(master_phil, local_overrides=override)

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    ###########################
    # Parameterise the models #
    ###########################

    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisationOrientation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)

    # Fix beam to the X-Z plane (imgCIF geometry)
    s0_param.set_fixed([True, False])

    # Fix crystal parameters
    #xluc_param.set_fixed([True, True, True, True, True, True])

    ########################################################################
    # Link model parameterisations together into a parameterisation of the #
    # prediction equation                                                  #
    ########################################################################

    pred_param = DetectorSpacePredictionParameterisation(
    mydetector, mybeam, mycrystal, mygonio, [det_param], [s0_param],
    [xlo_param], [xluc_param])

    print "The initial experimental geometry is:"
    print_model_geometry(mybeam, mydetector, mycrystal)
    print "\nInitial values of parameters are"
    msg = "Parameters: " + "%.5f " * len(pred_param)
    print msg % tuple(pred_param.get_param_vals())
    print

    return(mydetector, mygonio, mycrystal, mybeam,
           det_param, s0_param, xlo_param, xluc_param, pred_param)

def run(mydetector, mygonio, mycrystal, mybeam,
     det_param, s0_param, xlo_param, xluc_param,
     pred_param):

    #################################
    # Apply random parameter shifts #
    #################################

    # shift detector by normal deviate of sd 2.0 mm each translation and 4 mrad
    # each rotation
    det_p = det_param.get_param_vals()
    shift_det_p = random_param_shift(det_p, [2.0, 2.0, 2.0, 4.0, 4.0, 4.0])
    det_param.set_param_vals(shift_det_p)

    # rotate beam by normal deviate with sd 4 mrad. There is only one free axis!
    s0_p = s0_param.get_param_vals()
    shift_s0_p = random_param_shift(s0_p, [4.0])
    s0_param.set_param_vals(shift_s0_p)

    # rotate crystal by normal deviates with sd 4 mrad for each rotation.
    xlo_p = xlo_param.get_param_vals()
    shift_xlo_p = random_param_shift(xlo_p, [4.0, 4.0, 4.0])
    xlo_param.set_param_vals(shift_xlo_p)

    # change unit cell a bit (by normal deviates of 0.5 Angstrom length
    # upsets, 0.5 degree of gamma angle only)
    xluc_p_vals = xluc_param.get_param_vals()
    cell_params = mycrystal.get_unit_cell().parameters()
    shift_cell = random_param_shift(cell_params, [0.5, 0.5, 0.5,
                                                  0.0, 0.0, 0.5])
    new_uc = unit_cell(shift_cell)
    newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
    S = symmetrize_reduce_enlarge(mycrystal.get_space_group())
    S.set_orientation(orientation=newB)
    X = tuple([e * 1.e5 for e in S.forward_independent_parameters()])
    xluc_param.set_param_vals(X)

    target_param_values = tuple(pred_param.get_param_vals())

    #############################
    # Generate some reflections #
    #############################

    # All indices in a 2.0 Angstrom sphere
    resolution = 2.0
    index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                    space_group(space_group_symbols(1).hall()).type(), resolution)
    indices = index_generator.to_array()

    # Build a mock scan for a 180 degree sweep
    sf = scan_factory()
    myscan = sf.make_scan(image_range = (1,1800),
                          exposure_time = 0.1,
                          oscillation = (0, 0.1),
                          epochs = range(1800),
                          deg = True)
    sweep_range = myscan.get_oscillation_range(deg=False)
    temp = myscan.get_oscillation(deg=False)
    im_width = temp[1] - temp[0]
    assert sweep_range == (0., pi)
    assert approx_equal(im_width, 0.1 * pi / 180.)

    ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)
    obs_refs = ref_predictor.predict(indices)

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.
    px_size = mydetector[0].get_pixel_size()
    var_x = (px_size[0] / 2.)**2
    var_y = (px_size[1] / 2.)**2
    var_phi = (im_width / 2.)**2

    for ref in obs_refs:

        # calc and set the impact position, assuming all reflections
        # intersect panel 0.
        impacts = mydetector[0].get_ray_intersection(ref.beam_vector)
        ref.image_coord_mm = impacts

        # set the 'observed' centroids
        ref.centroid_position = ref.image_coord_mm + (ref.rotation_angle, )

        # set the centroid variance
        ref.centroid_variance = (var_x, var_y ,var_phi)

        # set the frame number, calculated from rotation angle
        ref.frame_number = myscan.get_image_index_from_angle(
            ref.rotation_angle, deg=False)

    ###############################
    # Undo known parameter shifts #
    ###############################

    s0_param.set_param_vals(s0_p)
    det_param.set_param_vals(det_p)
    xlo_param.set_param_vals(xlo_p)

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(obs_refs, mybeam, mygonio, sweep_range)

    ##############################
    # Set up the target function #
    ##############################

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(ref_predictor,
        mydetector, refman, pred_param, im_width)

    ################################
    # Set up the refinement engine #
    ################################

    #refiner = SimpleLBFGS(mytarget, pred_param)
    #refiner = LBFGScurvs(mytarget, pred_param)
    refiner = GaussNewtonIterations(mytarget, pred_param)
    refiner.run()

    msg = "%d %s " + "%.5f " * len(pred_param)
    subs = ((refiner.get_num_steps(), str(refiner._target_achieved)) +
            target_param_values)
    print msg % subs

    ##########################
    # Reset parameter values #
    ##########################

    s0_param.set_param_vals(s0_p)
    det_param.set_param_vals(det_p)
    xlo_param.set_param_vals(xlo_p)

    return refiner

if __name__ == "__main__":
    args = sys.argv[1:]
    try:
        seed = int(args[0])
    except IndexError:
        seed = 1
    if len(args) != 1: print "Usage:",sys.argv[0],"seed"

    (mydetector, mygonio, mycrystal, mybeam,
        det_param, s0_param, xlo_param, xluc_param,
        pred_param) = setup_models(seed)

    header = "Nsteps Completed " + "Param_%02d " * len(pred_param)
    print header % tuple(range(1, len(pred_param) + 1))

    for i in xrange(50):
        sys.stdout.flush()
        output = run(mydetector, mygonio, mycrystal, mybeam,
            det_param, s0_param, xlo_param, xluc_param,
            pred_param)
