#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Test analytical calculation of gradients of the target function versus finite
difference calculations"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
import random
from scitbx import matrix
from libtbx.phil import parse

# Experimental model builder
from setup_geometry import extract

# Model parameterisations
from dials.scratch.dgw.refinement.detector_parameters import \
    detector_parameterisation_single_sensor
from dials.scratch.dgw.refinement.source_parameters import \
    beam_parameterisation_orientation
from dials.scratch.dgw.refinement.crystal_parameters import \
    crystal_orientation_parameterisation, crystal_unit_cell_parameterisation
from dials.scratch.dgw.refinement import random_param_shift

# Reflection prediction
from dials.scratch.dgw.prediction import angle_predictor, impact_predictor
from rstbx.diffraction import full_sphere_indices
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.scratch.dgw.refinement.prediction_parameters import \
    detector_space_prediction_parameterisation

# Imports for the target function
from dials.scratch.dgw.refinement.target import \
    least_squares_positional_residual_with_rmsd_cutoff, reflection_manager

# Local functions
def random_direction_close_to(vector, sd = 0.5):
    return vector.rotate_around_origin(matrix.col(
                (random.random(),
                 random.random(),
                 random.random())).normalize(),
                 random.gauss(0, sd),  deg = True)

#############################
# Setup experimental models #
#############################

args = sys.argv[1:]
# make a small cell to speed up calculations
overrides = """geometry.parameters.crystal.a.length.range = 10 15
geometry.parameters.crystal.b.length.range = 10 15
geometry.parameters.crystal.c.length.range = 10 15"""

master_phil = parse("""
    include file geometry.params
    """, process_includes=True)

models = extract(master_phil, overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

print "Initial experimental model"
print "=========================="
print "beam s0       = (%.4f, %.4f, %.4f)" % mybeam.get_s0()
print "sensor origin = (%.4f, %.4f, %.4f)" % mydetector.sensors()[0].origin
print "sensor dir1   = (%.4f, %.4f, %.4f)" % mydetector.sensors()[0].dir1
print "sensor dir2   = (%.4f, %.4f, %.4f)" % mydetector.sensors()[0].dir2
uc = mycrystal.get_unit_cell()
print "crystal unit cell = %.4f, %.4f, %.4f, %.4f, %.4f, %.4f" % uc.parameters()
print "crystal orientation matrix U ="
print mycrystal.get_U().round(4)
print

###########################
# Parameterise the models #
###########################

det_param = detector_parameterisation_single_sensor(mydetector.sensors()[0])
s0_param = beam_parameterisation_orientation(mybeam)
xlo_param = crystal_orientation_parameterisation(mycrystal)
xluc_param = crystal_unit_cell_parameterisation(mycrystal) # dummy, does nothing

########################################################################
# Link model parameterisations together into a parameterisation of the #
# prediction equation                                                  #
########################################################################

# Testing the new 'coupled' version here
pred_param = detector_space_prediction_parameterisation(
    mydetector, mybeam, mycrystal, mygonio, [det_param], [s0_param],
    [xlo_param], [xluc_param])

################################
# Apply known parameter shifts #
################################

# shift detector by 0.2 mm each translation and 2 mrad each rotation
det_p_vals = det_param.get_p()
p_vals = [a + b for a, b in zip(det_p_vals,
                                 [2.0, 2.0, 2.0, 2.0, 2.0, 2.0])]
det_param.set_p(p_vals)

# shift beam by 2 mrad in one axis
s0_p_vals = s0_param.get_p()
p_vals = list(s0_p_vals)
p_vals[1] += 2.0
s0_param.set_p(p_vals)

# rotate crystal a bit (=4 mrad each rotation)
#p_vals = xlo_param.get_p()
#p_vals2 = random_param_shift(p_vals, [0.004, 0.004, 0.004])
#p_vals2 = [0.004, 0.004, 0.004]
#xlo_param.set_p(p_vals2)
xlo_p_vals = xlo_param.get_p()
p_vals = [a + b for a, b in zip(xlo_p_vals, [2.0, 2.0, 2.0])]
xlo_param.set_p(p_vals)

print "Offsetting initial model"
print "========================"
print "beam s0       = (%.4f, %.4f, %.4f)" % mybeam.get_s0()
print "sensor origin = (%.4f, %.4f, %.4f)" % mydetector.sensors()[0].origin
print "sensor dir1   = (%.4f, %.4f, %.4f)" % mydetector.sensors()[0].dir1
print "sensor dir2   = (%.4f, %.4f, %.4f)" % mydetector.sensors()[0].dir2
uc = mycrystal.get_unit_cell()
print "crystal unit cell = %.4f, %.4f, %.4f, %.4f, %.4f, %.4f" % uc.parameters()
print "crystal orientation matrix U ="
print mycrystal.get_U().round(4)
print

#############################
# Generate some reflections #
#############################

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
indices = full_sphere_indices(
    unit_cell = mycrystal.get_unit_cell(),
    resolution_limit = resolution,
    space_group = space_group(space_group_symbols(1).hall()))

# Select those that are excited in a 180 degree sweep and get their angles
UB = mycrystal.get_U() * mycrystal.get_B()
ap = angle_predictor(mycrystal, mybeam, mygonio, resolution)
obs_indices, obs_angles = ap.observed_indices_and_angles_from_angle_range(
    phi_start_rad = 0.0, phi_end_rad = pi, indices = indices)
print "Generating reflections"
print "======================"
print "Total number of reflections excited", len(obs_indices)

# Project positions on camera
ip = impact_predictor(mydetector, mygonio, mybeam, mycrystal)
hkls, d1s, d2s, angles, svecs = ip.predict(obs_indices.as_vec3_double(),
                                       obs_angles)
print "Total number of observations made", len(hkls), "\n"

# Invent some uncertainties
im_width = 0.1 * pi / 180.
sigd1s = [mydetector.px_size_fast() / 2.] * len(hkls)
sigd2s = [mydetector.px_size_slow() / 2.] * len(hkls)
sigangles = [im_width / 2.] * len(hkls)

###############################
# Undo known parameter shifts #
###############################

s0_param.set_p(s0_p_vals)
det_param.set_p(det_p_vals)
xlo_param.set_p(xlo_p_vals)

print "Resetting to initial model"
print "=========================="
print "Initial values of parameters are"
msg = "Parameters: " + "%.5f " * len(pred_param)
print msg % tuple(pred_param.get_p()), "\n"

#####################################
# Select reflections for refinement #
#####################################

refman = reflection_manager(hkls, svecs,
                        d1s, sigd1s,
                        d2s, sigd2s,
                        angles, sigangles,
                        mybeam, mygonio)

##############################
# Set up the target function #
##############################

mytarget = least_squares_positional_residual_with_rmsd_cutoff(
    refman, ap, ip, pred_param, mydetector.px_size_fast(),
    mydetector.px_size_slow(), im_width)

# get the functional and gradients
mytarget.predict()
L, dL_dp = mytarget.compute_functional_and_gradients()
try:
    curvs = mytarget.curvatures()
except AttributeError:
    curvs = None

print "Calculated gradients"
print "===================="
print "Target L           = %.6f" % L
msg = "Gradients dL/dp    = " + "%.6f " * len(dL_dp)
print msg % tuple(dL_dp)
if curvs:
    msg = "Curvatures d2L/dp2 = " + "%.6f " * len(curvs)
    print msg % tuple(curvs)
print

####################################
# Do FD calculation for comparison #
####################################

# function for calculating finite difference gradients of the target function
def get_fd_gradients(target, pred_param, deltas):
    '''Calculate centered finite difference gradients for each of the
    parameters of the target function.

    "deltas" must be a sequence of the same length as the parameter list, and
    contains the step size for the difference calculations for each parameter.
    '''

    p_vals = pred_param.get_p()
    assert len(deltas) == len(p_vals)
    fd_grad = []
    fd_curvs = []

    try:
        do_curvs = callable(target.curvatures)
    except AttributeError:
        do_curvs = False

    print "Parameter",
    for i in range(len(deltas)):
        val = p_vals[i]

        print i,
        p_vals[i] -= deltas[i] / 2.
        pred_param.set_p(p_vals)
        target.predict()

        rev_state = target.compute_functional_and_gradients()

        p_vals[i] += deltas[i]
        pred_param.set_p(p_vals)

        target.predict()

        fwd_state = target.compute_functional_and_gradients()

        # finite difference estimation of first derivatives
        fd_grad.append((fwd_state[0] - rev_state[0]) / deltas[i])

        # finite difference estimation of curvatures, using the analytical
        # first derivatives
        if do_curvs:
            fd_curvs.append((fwd_state[1][i] - rev_state[1][i]) / deltas[i])

        # set parameter back to centred value
        p_vals[i] = val

    print "\n"

    # return to the initial state
    pred_param.set_p(p_vals)

    return fd_grad, fd_curvs

print "Finite difference gradients"
print "==========================="
fdgrads = get_fd_gradients(mytarget, pred_param, [1.e-7] * len(pred_param))
msg = "FD gradients dL[fd]/dp          = " + "%.6f " * len(fdgrads[0])
print msg % tuple(fdgrads[0])
diffs = [a - b for a, b in zip(dL_dp, fdgrads[0])]
msg = "dL/dp - dL[fd]/dp               = " + "%.6f " * len(diffs)
print msg % tuple(diffs)
print "Normalised differences:"
msg = "(dL/dp - dL[fd]/dp) / dL[fd]/dp = " + "%.6f " * len(diffs)
print msg % tuple([a / b for a, b in zip(diffs, fdgrads[0])])
print

if curvs:
    print "Finite difference curvatures"
    print "============================"
    msg = "FD curvatures d2L[fd]/dp2             = " + "%.6f " * len(fdgrads[1])
    print msg % tuple(fdgrads[1])
    diffs = [a - b for a, b in zip(curvs, fdgrads[1])]
    msg = "d2L/dp2 - d2L[fd]/dp2                 = " + "%.6f " * len(diffs)
    print msg % tuple(diffs)
    print "Normalised differences:"
    msg = "(d2L/dp2 - d2L[fd]/dp2) / d2L[fd]/dp2 = " + "%.6f " * len(diffs)
    print msg % tuple([a / b for a, b in zip(diffs, fdgrads[1])])
