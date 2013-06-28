"""Auxiliary functions for the refinement package"""

from __future__ import division
from math import sin, cos, sqrt
from scitbx import matrix
import random

# flex is required before the boost python import to avoid an import error.
# Ignore complaints by libtbx.find_clutter
from cctbx.array_family import flex
from dials_refinement_helpers_ext import *


def dR_from_axis_and_angle(axis, angle, deg=False):
    '''return the first derivative of a rotation matrix specified by its
    axis and angle'''

    # NB it is inefficient to do this separately from the calculation of
    # the rotation matrix itself, but it seems the Python interface to
    # scitbx does not have a suitable function. It might perhaps be
    # useful to write one, which could come straight from David Thomas'
    # RTMATS (present in Mosflm and MADNES).

    # NB RTMATS does calculation for a clockwise rotation of a vector
    # whereas axis_and_angle_as_r3_rotation_matrix does anticlockwise
    # rotation. Therefore flip the axis in here compared with
    # RTMATS in order to match the axis_and_angle_as_r3_rotation_matrix
    # convention

    assert axis.n in ((3,1), (1,3))
    if (deg): angle *= pi/180
    axis *= -1.
    ca, sa  = cos(angle), sin(angle)

    return(matrix.sqr((sa * axis[0] * axis[0] - sa ,
                    sa * axis[0] * axis[1] + ca * axis[2],
                    sa * axis[0] * axis[2] - ca * axis[1],
                    sa * axis[1] * axis[0] - ca * axis[2],
                    sa * axis[1] * axis[1] - sa,
                    sa * axis[1] * axis[2] + ca * axis[0],
                    sa * axis[2] * axis[0] + ca * axis[1],
                    sa * axis[2] * axis[1] - ca * axis[0],
                    sa * axis[2] * axis[2] - sa)))

def random_param_shift(vals, sigmas):
    '''Add a random (normal) shift to a parameter set, for testing'''

    assert len(vals) == len(sigmas)
    shifts = [random.gauss(0, sd) for sd in sigmas]
    newvals = [(x + y) for x, y in zip(vals, shifts)]

    return newvals

def get_fd_gradients(mp, deltas):
    '''Calculate centered finite difference gradients for each of the
    parameters of the model parameterisation mp.

    "deltas" must be a sequence of the same length as the parameter list, and
    contains the step size for the difference calculations for each parameter.
    '''

    #state = matrix.sqr(mp._models[0].d)

    p_vals = mp.get_p()
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):

        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.
        mp.set_p(p_vals)
        rev_state = mp.get_state()

        p_vals[i] += deltas[i]
        mp.set_p(p_vals)
        fwd_state = mp.get_state()

        fd_grad.append((fwd_state - rev_state) / deltas[i])

        p_vals[i] = val

    # return to the initial state
    mp.set_p(p_vals)

    return fd_grad

def print_model_geometry(beam = None, detector = None, crystal = None):

    if beam:
        print "beam s0 = (%.4f, %.4f, %.4f)" % beam.get_s0()
    if detector:
        print "sensor origin = (%.4f, %.4f, %.4f)" % detector[0].get_origin()
        print "sensor dir1 = (%.4f, %.4f, %.4f)" % detector[0].get_fast_axis()
        print "sensor dir2 = (%.4f, %.4f, %.4f)" % detector[0].get_slow_axis()
    if crystal:
        uc = crystal.get_unit_cell()
        print "crystal unit cell = %.4f, %.4f, %.4f, %.4f, %.4f, %.4f" % uc.parameters()
        print "crystal orientation matrix U ="
        print crystal.get_U().round(4)

def print_grads(grad_list):
    for i, grad in enumerate(grad_list):
        print ("Param %02d. Gradients: "
               "%.5f, %.5f, %.5f" % ((i,) + tuple(grad)))

def refine(beam, goniometer, crystal, detector, scan,
           reflections, verbosity = 0, fix_cell = False,
           scan_varying=False):

    """Simple refinement interface for the centroid refinement sprint"""

    # Reflection prediction
    from dials.algorithms.refinement.prediction import ReflectionPredictor

    # Model parameterisations
    from dials.algorithms.refinement.parameterisation.detector_parameters import \
        DetectorParameterisationSinglePanel
    from dials.algorithms.refinement.parameterisation.source_parameters import \
        BeamParameterisationOrientation
    from dials.algorithms.refinement.parameterisation.crystal_parameters import \
        CrystalOrientationParameterisation, CrystalUnitCellParameterisation

    # Symmetry constrained parameterisation for the unit cell
    #from cctbx.uctbx import unit_cell
    #from rstbx.symmetry.constraints.parameter_reduction import \
    #    symmetrize_reduce_enlarge

    # Parameterisation of the prediction equation
    from dials.algorithms.refinement.parameterisation.prediction_parameters import \
        DetectorSpacePredictionParameterisation

    # Imports for the target function
    from dials.algorithms.refinement.target import \
        LeastSquaresPositionalResidualWithRmsdCutoff, ReflectionManager

    # Import the refinement engine
    from dials.algorithms.refinement.engine import GaussNewtonIterations

    # pull out data needed for refinement
    temp = [(ref.miller_index, ref.entering, ref.frame_number,
             ref.rotation_angle, matrix.col(ref.beam_vector),
             ref.image_coord_mm, ref.centroid_variance) \
                for ref in reflections]
    hkls, enterings, frames, angles, svecs, intersects, variances = zip(*temp)

    # tease apart tuples to separate lists
    d1s, d2s = zip(*intersects)
    var_d1s, var_d2s, var_angles = zip(*variances)

    # change variances to sigmas
    sig_d1s = [sqrt(e) for e in var_d1s]
    sig_d2s = [sqrt(e) for e in var_d2s]
    sig_angles = [sqrt(e) for e in var_angles]

    assert len(hkls) == len(svecs) == len(d1s) == len(d2s) == \
           len(sig_d2s) == len(angles) == len(sig_angles)

    image_width = scan.get_oscillation(deg=False)[1]
    sweep_range = scan.get_oscillation_range(deg=False)
    ref_predictor = ReflectionPredictor(crystal, beam, goniometer, sweep_range)

    ###########################
    # Parameterise the models #
    ###########################

    det_param = DetectorParameterisationSinglePanel(detector)
    s0_param = BeamParameterisationOrientation(beam, goniometer)
    xlo_param = CrystalOrientationParameterisation(crystal)
    xluc_param = CrystalUnitCellParameterisation(crystal)

    # Fix beam to the X-Z plane (imgCIF geometry)
    s0_param.set_fixed([True, False])

    # Fix cell if requested
    if fix_cell:
        xluc_param.set_fixed([True] * xluc_param.num_free())

    ########################################################################
    # Link model parameterisations together into a parameterisation of the #
    # prediction equation                                                  #
    ########################################################################

    pred_param = DetectorSpacePredictionParameterisation(
    detector, beam, crystal, goniometer, [det_param], [s0_param],
    [xlo_param], [xluc_param])

    if verbosity > 1:
        print "Prediction equation parameterisation built\n"
        print "Parameter order : name mapping"
        for i, e in enumerate(pred_param.get_p_names()):
            print "Parameter %03d : " % i + e
        print

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(hkls, enterings, frames, svecs,
                               d1s, sig_d1s,
                               d2s, sig_d2s,
                               angles, sig_angles,
                               beam, goniometer, scan, verbosity)

    if verbosity > 1: print "Reflection manager built\n"

    ##############################
    # Set up the target function #
    ##############################

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
        ref_predictor, detector, refman, pred_param, image_width)

    if verbosity > 1: print "Target function built\n"

    ################################
    # Set up the refinement engine #
    ################################

    refiner = GaussNewtonIterations(mytarget, pred_param, log=None,
                                    verbosity=verbosity)

    if verbosity > 1: print "Refinement engine built\n"

    ###################################
    # Do refinement and return models #
    ###################################

    refiner.run()

    #############################################################
    # Do a second round of refinement with scan-varying models? #
    #############################################################

    if scan_varying: scan_varying_refine(beam, goniometer, crystal,
           detector, image_width, scan,
           hkls, enterings, frames,
           svecs, d1s, sig_d1s, d2s, sig_d2s, angles, sig_angles,
           verbosity, fix_cell = fix_cell)

    # Return models (crystal not updated with scan-varying U, B)
    return(beam, goniometer, crystal, detector)


def scan_varying_refine(beam, goniometer, crystal, detector,
           image_width, scan,
           hkls,  enterings, frames,
           svecs, d1s, sigd1s, d2s, sigd2s, angles, sigangles,
           verbosity = 0, fix_cell = False):

    """experimental refinement function for scan-varying refinement"""

    # Reflection prediction
    from dials.algorithms.refinement.prediction import ReflectionPredictor

    # Model parameterisations
    from dials.algorithms.refinement.parameterisation.detector_parameters import \
        DetectorParameterisationSinglePanel
    from dials.algorithms.refinement.parameterisation.source_parameters import \
        BeamParameterisationOrientation
    from dials.algorithms.refinement.parameterisation.\
        scan_varying_crystal_parameters import \
            ScanVaryingCrystalOrientationParameterisation, \
            ScanVaryingCrystalUnitCellParameterisation

    # Symmetry constrained parameterisation for the unit cell
    #from cctbx.uctbx import unit_cell
    #from rstbx.symmetry.constraints.parameter_reduction import \
    #    symmetrize_reduce_enlarge

    # Parameterisation of the prediction equation
    from dials.algorithms.refinement.parameterisation.\
        scan_varying_prediction_parameters \
            import VaryingCrystalPredictionParameterisation

    # Imports for the target function
    from dials.scratch.dgw.refinement.scan_varying_target import \
        LeastSquaresPositionalResidualWithRmsdCutoff
    from dials.scratch.dgw.refinement.scan_varying_target import \
        ReflectionManager

    # Import the refinement engine
    from dials.algorithms.refinement.engine import GaussNewtonIterations

    sweep_range = scan.get_oscillation_range(deg=False)
    ref_predictor = ReflectionPredictor(crystal, beam, goniometer, sweep_range)

    ###########################
    # Parameterise the models #
    ###########################

    det_param = DetectorParameterisationSinglePanel(detector)
    s0_param = BeamParameterisationOrientation(beam, goniometer)
    xlo_param = ScanVaryingCrystalOrientationParameterisation(
            crystal, scan.get_image_range(), 5)
    xluc_param = ScanVaryingCrystalUnitCellParameterisation(
            crystal, scan.get_image_range(), 5)

    # Fix beam to the X-Z plane (imgCIF geometry)
    s0_param.set_fixed([True, False])

    # Fix cell if requested
    if fix_cell:
        xluc_param.set_fixed([True] * xluc_param.num_free())

    ########################################################################
    # Link model parameterisations together into a parameterisation of the #
    # prediction equation                                                  #
    ########################################################################

    pred_param = VaryingCrystalPredictionParameterisation(
        detector, beam, crystal, goniometer, [det_param], [s0_param],
        [xlo_param], [xluc_param])

    if verbosity > 1:
        print "Prediction equation parameterisation built\n"
        print "Parameter order:name mapping"
        for i, e in enumerate(pred_param.get_p_names()):
            print "Parameter %03d : " % i + e
        print

    #####################################
    # Select reflections for refinement #
    #####################################

    refman = ReflectionManager(hkls, enterings, frames,
                            svecs,
                            d1s, sigd1s,
                            d2s, sigd2s,
                            angles, sigangles,
                            beam, goniometer, scan, verbosity,
                            nref_per_degree = 50)

    if verbosity > 1:
        print "Reflection manager built\n"
        print "Working set size = %d observations" % refman.get_sample_size()

    ##############################
    # Set up the target function #
    ##############################

    mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
        ref_predictor, detector, refman, pred_param, image_width)

    if verbosity > 1: print "Target function built\n"

    ################################
    # Set up the refinement engine #
    ################################

    refiner = GaussNewtonIterations(mytarget, pred_param, log=None,
                                    verbosity=verbosity)

    if verbosity > 1: print "Refinement engine built\n"

    ###################################
    # Do refinement and return models #
    ###################################

    refiner.run()

    # NB scan-independent models set by side-effect, scan-varying
    # crystal model not yet set at all
    return
