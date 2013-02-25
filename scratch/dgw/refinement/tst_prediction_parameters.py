#!/usr/bin/env python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

from __future__ import division

#### Python and general cctbx imports

import sys
import random
from math import pi
from scitbx import matrix
#from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

##### Import model builder
from setup_geometry import extract

#### Import model parameterisations

from dials.scratch.dgw.refinement.prediction_parameters import \
    detector_space_prediction_parameterisation
from dials.scratch.dgw.refinement.detector_parameters import \
    detector_parameterisation_single_sensor
from dials.scratch.dgw.refinement.source_parameters import \
    source_parameterisation_orientation
from dials.scratch.dgw.refinement.crystal_parameters import \
    crystal_orientation_parameterisation, crystal_unit_cell_parameterisation

#### Local functions

def print_grads(grad_list):
    for i, grad in enumerate(grad_list):
        print ("Param %02d. Gradients: "
               "%.5f, %.5f, %.5f" % ((i,) + tuple(grad)))

# Functions required for finite difference calculations

def get_state(gon, src, xl, det, hkl, angle, angle_predictor):
    '''reflection prediction for the current state of the models'''

    # predict the rotation angles
    obs_ang = angle_predictor.predict(hkl)

    # select which is nearest the observed angle
    deltas = [abs(x - angle) for x in obs_ang]
    new_angle = obs_ang[deltas.index(min(deltas))]

    rp = reflection_prediction(gon.get_axis(),
                               src.get_s0(),
                               xl.get_U() * xl.get_B(),
                               det.sensors()[0])

    # cache prediction for one hkl at angle
    if rp(hkl, new_angle):
        return matrix.col(rp.get_prediction() + (new_angle, ))
    else: return None

def get_fd_gradients(pred_param, hkl, phi, angle_predictor, deltas):
    '''Calculate centered finite difference gradients for each of the
    parameters of the prediction parameterisation object, for reflection hkl
    at angle phi.

    "deltas" must be a sequence of the same length as the parameter list,
    and contains the step size for the difference calculations for each
    parameter.'''

    gon = pred_param._gonio
    src = pred_param._source
    det = pred_param._detector
    xl = pred_param._crystal
    ap = angle_predictor

    p_vals = pred_param.get_p()
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):

        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.
        pred_param.set_p(p_vals)
        rev_state = get_state(gon, src, xl, det, hkl, phi, ap)

        p_vals[i] += deltas[i]
        pred_param.set_p(p_vals)
        fwd_state = get_state(gon, src, xl, det, hkl, phi, ap)

        fd_grad.append((fwd_state - rev_state) / deltas[i])
        p_vals[i] = val

    # return to the initial state
    pred_param.set_p(p_vals)

    return fd_grad

from time import time

start_time = time()

#### Create models

args = sys.argv[1:]
overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

models = extract(overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mysource = models.source

#### Create parameterisations of these models

det_param = detector_parameterisation_single_sensor(mydetector.sensors()[0])
src_param = source_parameterisation_orientation(mysource)
xlo_param = crystal_orientation_parameterisation(mycrystal)
xluc_param = crystal_unit_cell_parameterisation(mycrystal)

#### Unit tests

# Build a prediction parameterisation with a single detector model
pred_param = detector_space_prediction_parameterisation(mydetector,
             mysource, mycrystal, mygonio, [det_param])

# Check the accessors
assert len(pred_param) == 6
for (a, b) in zip(pred_param.get_p(), det_param.get_p()):
    assert a==b

pred_param.set_p([100., 1.0, 1.0, 0., 0., 0.])
for (a, b) in zip(pred_param.get_p(), det_param.get_p()):
    assert a==b

# Build a full global parameterisation with detector and source
# parameterised and not-fully-functional crystal parameterisations
pred_param = detector_space_prediction_parameterisation(
    mydetector,  mysource, mycrystal, mygonio, [det_param], [src_param],
    [xlo_param], [xluc_param])

# Generate some indices
from dials.scratch.dgw.prediction import angle_predictor
from rstbx.diffraction import reflection_prediction
from rstbx.diffraction import full_sphere_indices
from cctbx.sgtbx import space_group, space_group_symbols
resolution = 2.0
indices = full_sphere_indices(
    unit_cell = mycrystal.get_unit_cell(),
    resolution_limit = resolution,
    space_group = space_group(space_group_symbols(1).hall()))

# Generate list of phi values
UB = mycrystal.get_U() * mycrystal.get_B()
ap = angle_predictor(mycrystal, mysource, mygonio, resolution)
obs_indices, obs_angles = ap.observed_indices_and_angles_from_angle_range(
    phi_start_rad = 0.0, phi_end_rad = pi/5., indices = indices)

# Project positions on camera
rp = reflection_prediction(mygonio.get_axis(), mysource.get_s0(), UB,
                           mydetector.sensors()[0])
hkls, d1s, d2s, angles, s_dirs = rp.predict(obs_indices.as_vec3_double(),
                                       obs_angles)

# Test get_state for the first reflection
tmp = get_state(mygonio, mysource, mycrystal, mydetector, hkls[0],
                angles[0], ap)
for (a, b) in zip(tmp, (d1s[0], d2s[0], angles[0])):
    assert a == b

# Compare analytical and finite difference gradients for up to 50 randomly
# selected reflections.
# NB, reflections that just touch the Ewald sphere have large derivatives of phi
# wrt some parameters (approching infinity). Such reflections are likely to fail
# the comparison against FD gradients. Andrew Leslie points out that these
# reflections are not even integrated, because of susceptibility to errors in
# the Lorentz correction. To avoid failures of this test, we exclude reflections
# that are closer than the largest reciprocal lattice vector length from the
# plane of the rotation axis and beam direction.
selection = random.sample(xrange(len(hkls)), min(len(hkls), 50))
uc = mycrystal.get_unit_cell()
exclusion_limit = max(uc.reciprocal_parameters()[0:3])

verbose = False
for iref in selection:
    hkl, s_dir, angle = hkls[iref], s_dirs[iref], angles[iref]

    # Beware! s_dirs[0] has been normalised. I need the proper length vector s
    s = matrix.col(s_dir) / mysource.get_wavelength()

    # get analytical gradients
    an_grads = pred_param.get_gradients(hkl, s, angle)

    # Reflections close to being in the plane of the rotation axis and beam
    # give very large gradients of phi, and may fail the test versus finite
    # difference gradients. Detect and exclude very large gradients
    s0 = mysource.get_s0()
    r = s - s0
    e_s0_plane_norm = mygonio.get_axis().cross(s0).normalize()
    r_dist_from_plane = abs(r.dot(e_s0_plane_norm))
    if r_dist_from_plane <= exclusion_limit:
        continue

    fd_grads = get_fd_gradients(pred_param, hkl, angle, ap,
                                [1.e-7] * len(pred_param))

    if verbose:
        print hkl
        for g in an_grads: print g
        for g in fd_grads: print g
    for i, (an_grad, fd_grad) in enumerate(zip(an_grads, fd_grads)):
        for a, b in zip(an_grad, fd_grad):
            try:
                if abs(b) > 10.:
                    assert approx_equal((a - b) / b, 0., eps = 5.e-6)
                else:
                    assert approx_equal(a - b, 0., eps = 5.e-6)
            except AssertionError:
                print "Failure for parameter number %d" %i
                print "Analytical derivatives: %.6f, %.6f, %.6f" % tuple(an_grad)
                print "Finite derivatives: %.6f, %.6f, %.6f" % tuple(fd_grad)
                finish_time = time()
                print "Time Taken: ",finish_time - start_time
                raise


# if we got this far,
print "OK"

