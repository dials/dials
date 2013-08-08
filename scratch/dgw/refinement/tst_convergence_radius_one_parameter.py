#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Investigate the convergence radius of refinement by multiple refinement runs
changing one parameter at a time"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from scitbx import matrix
from libtbx.phil import parse

# Get class to build experimental models
from setup_geometry import Extract

# Model parameterisations
from dials.scratch.dgw.refinement.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.scratch.dgw.refinement.beam_parameters import \
    BeamParameterisationOrientation
from dials.scratch.dgw.refinement.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation
from dials.scratch.dgw.refinement import random_param_shift

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.scratch.dgw.prediction import ReflectionPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.scratch.dgw.refinement.prediction_parameters import \
    DetectorSpacePredictionParameterisation

# Imports for the target function
from dials.scratch.dgw.refinement.target_old import \
    LeastSquaresPositionalResidualWithRmsdCutoff, ReflectionManager

# Import the refinement engine
from dials.scratch.dgw.refinement.engine import SimpleLBFGS, LBFGScurvs

# Import helper functions
from dials.scratch.dgw.refinement import print_model_geometry

def setup_models(seed):

    #############################
    # Setup experimental models #
    #############################

    override = "geometry.parameters.random_seed=" + str(seed)
    print override

    master_phil = parse("""
    include file geometry.params
    include file minimiser.params
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
    xluc_param.set_fixed([True, True, True, True, True, True])

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
    print msg % tuple(pred_param.get_p())
    print

    return(mydetector, mygonio, mycrystal, mybeam,
           det_param, s0_param, xlo_param, xluc_param, pred_param)

def run(mydetector, mygonio, mycrystal, mybeam,
     det_param, s0_param, xlo_param, xluc_param,
     pred_param):

    #################################
    # Apply random parameter shifts #
    #################################

    # shift one detector param by normal deviate of sd 10.0 mm translation
    det_p = det_param.get_p()
    shift_det_p = random_param_shift(det_p, [10.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    det_param.set_p(shift_det_p)

    # rotate beam by normal deviate with sd 4 mrad. There is only one free axis!
    s0_p = s0_param.get_p()
    #shift_s0_p = random_param_shift(s0_p, [4.0])
    #s0_param.set_p(shift_s0_p)

    # rotate crystal by normal deviates with sd 4 mrad for each rotation.
    xlo_p = xlo_param.get_p()
    #shift_xlo_p = random_param_shift(xlo_p, [4.0, 4.0, 4.0])
    #xlo_param.set_p(shift_xlo_p)

    target_param_values = tuple(pred_param.get_p())

    #############################
    # Generate some reflections #
    #############################

    # All indices in a 2.0 Angstrom sphere
    resolution = 2.0
    index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                    space_group(space_group_symbols(1).hall()).type(), resolution)
    indices = index_generator.to_array()

    # Select those that are excited in a 180 degree sweep and get their angles
    UB = mycrystal.get_U() * mycrystal.get_B()
    sweep_range = (0., pi)
    ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)
    obs_refs = ref_predictor.predict(indices)

    # Pull out reflection data as lists
    temp = [(ref.miller_index, ref.rotation_angle,
             matrix.col(ref.beam_vector)) for ref in obs_refs]
    hkls, angles, svecs = zip(*temp)

    # Project positions on camera
    # currently assume all reflections intersect panel 0
    impacts = [mydetector[0].get_ray_intersection(
                            ref.beam_vector) for ref in obs_refs]
    d1s, d2s = zip(*impacts)

    # Invent some uncertainties
    im_width = 0.1 * pi / 180.
    px_size = mydetector.get_pixel_size()
    sigd1s = [px_size[0] / 2.] * len(hkls)
    sigd2s = [px_size[1] / 2.] * len(hkls)
    sigangles = [im_width / 2.] * len(hkls)

    ###############################
    # Undo known parameter shifts #
    ###############################

    s0_param.set_p(s0_p)
    det_param.set_p(det_p)
    xlo_param.set_p(xlo_p)

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(hkls, svecs,
                            d1s, sigd1s,
                            d2s, sigd2s,
                            angles, sigangles,
                            mybeam, mygonio)

    ##############################
    # Set up the target function #
    ##############################

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
        refman, ref_predictor, mydetector, pred_param, im_width)

    ################################
    # Set up the refinement engine #
    ################################

    #refiner = SimpleLBFGS(mytarget, pred_param)
    refiner = LBFGScurvs(mytarget, pred_param)
    refiner.run()

    msg = "%d %s " + "%.5f " * len(pred_param)
    subs = ((refiner._step, str(refiner._target_achieved)) +
            target_param_values)
    print msg % subs

    ##########################
    # Reset parameter values #
    ##########################

    s0_param.set_p(s0_p)
    det_param.set_p(det_p)
    xlo_param.set_p(xlo_p)

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

    for i in xrange(100):
        sys.stdout.flush()
        output = run(mydetector, mygonio, mycrystal, mybeam,
            det_param, s0_param, xlo_param, xluc_param,
            pred_param)
