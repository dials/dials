#!/usr/bin/env python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division

#### Python and general cctbx imports

import sys
import random
from math import pi
from scitbx import matrix
from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.test_utils import approx_equal
from libtbx.phil import parse

##### Import model builders

from setup_geometry import Extract
from dxtbx.model.scan import scan_factory

##### Imports for reflection prediction

from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ReflectionPredictor

#### Import model parameterisations

from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters import \
    VaryingCrystalPredictionParameterisation
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters import \
    ScanVaryingCrystalOrientationParameterisation, \
    ScanVaryingCrystalUnitCellParameterisation

#### Import helper functions

from dials.algorithms.refinement import random_param_shift, print_grads

#### Local functions

# Functions required for finite difference calculations

def get_state(det, hkl, UB, angle, reflection_predictor):
    """reflection prediction for the current state of the models"""

    # update reflection_predictor with latest geometry
    reflection_predictor.update()

    # predict for this hkl
    refs = reflection_predictor.predict(hkl, UB)

    # select which is nearest the observed angle
    deltas = [abs(ref.rotation_angle - angle) for ref in refs]
    new_ref = refs[deltas.index(min(deltas))]

    # predict impact position
    impact = det[0].get_ray_intersection(new_ref.beam_vector)

    return matrix.col((impact) + (new_ref.rotation_angle,))

def get_fd_gradients(pred_param, hkl, phi, frame, reflection_predictor,
                     deltas):
    """Calculate centered finite difference gradients for each of the
    parameters of the prediction parameterisation object, for reflection
    hkl at angle phi.

    "deltas" must be a sequence of the same length as the parameter list,
    and contains the step size for the difference calculations for each
    parameter."""

    gon = pred_param._gonio
    src = pred_param._beam
    det = pred_param._detector
    xl = pred_param._crystal
    rp = reflection_predictor

    p_vals = pred_param.get_param_vals()
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):

        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.
        pred_param.set_param_vals(p_vals)

        # get UB for the current frame
        xlo_param.compose(frame)
        xluc_param.compose(frame)
        U = xlo_param.get_state()
        B = xluc_param.get_state()
        UB = U * B

        rev_state = get_state(det, hkl, UB, phi, rp)

        p_vals[i] += deltas[i]
        pred_param.set_param_vals(p_vals)

        # get UB for the current frame
        xlo_param.compose(frame)
        xluc_param.compose(frame)
        U = xlo_param.get_state()
        B = xluc_param.get_state()
        UB = U * B

        fwd_state = get_state(det, hkl, UB, phi, rp)

        fd_grad.append((fwd_state - rev_state) / deltas[i])
        p_vals[i] = val

    # return to the initial state
    pred_param.set_param_vals(p_vals)

    return fd_grad

from time import time

start_time = time()

#### Create models and parameterisations

args = sys.argv[1:]
overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    """, process_includes=True)

models = Extract(master_phil, overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

# Make a scan of 1-360 * 0.5 deg images
sf = scan_factory()
myscan = sf.make_scan((1,360), 0.5, (0, 0.5), range(360))
print myscan

# Create parameterisations of these models, with 5 samples for the
# scan-varying crystal parameterisations

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, mygonio)
xlo_param = ScanVaryingCrystalOrientationParameterisation(
        mycrystal, myscan.get_image_range(), 5)
xluc_param = ScanVaryingCrystalUnitCellParameterisation(
        mycrystal, myscan.get_image_range(), 5)

#### Cause the crystal U and B to vary over the scan

# Vary orientation angles by ~1.0 mrad each checkpoint
p_vals = xlo_param.get_param_vals()
sigmas = [1.0] * len(p_vals)
new_vals = random_param_shift(p_vals, sigmas)
xlo_param.set_param_vals(new_vals)

# Vary unit cell parameters, on order of 1% of the initial metrical
# matrix parameters
p_vals = xluc_param.get_param_vals()
sigmas = [0.01 * p for p in p_vals]
new_vals = random_param_shift(p_vals, sigmas)
xluc_param.set_param_vals(new_vals)

#### Unit tests

# Build a prediction equation parameterisation
pred_param = VaryingCrystalPredictionParameterisation(
    mydetector, mybeam, mycrystal, mygonio, [det_param], [s0_param],
    [xlo_param], [xluc_param])

# Generate some indices
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                      space_group(space_group_symbols(1).hall()).type(),
                      resolution)
indices = index_generator.to_array()

# Generate list of reflections
UB = mycrystal.get_U() * mycrystal.get_B()
ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio,
                                myscan.get_oscillation_range(deg=False))
ref_list = ref_predictor.predict(indices)

# Pull out lists of required reflection data
temp = [(ref.miller_index, ref.rotation_angle,
         matrix.col(ref.beam_vector)) for ref in ref_list]
hkls, angles, s_vecs = zip(*temp)

# convert angles to image number
frames = map(lambda x: myscan.get_image_index_from_angle(x, deg=False), angles)
print max(frames), min(frames), sum(frames) / len(frames)

# Project positions on camera. Currently assuming all reflections
# intersect panel 0
panel_id = 0
impacts = [mydetector[panel_id].get_ray_intersection(
                        ref.beam_vector) for ref in ref_list]
d1s, d2s = zip(*impacts)

# Test get_state for the first reflection
tmp = get_state(mydetector, hkls[0], UB, angles[0], ref_predictor)
for (a, b) in zip(tmp, (d1s[0], d2s[0], angles[0])):
    assert a == b

# Compare analytical and finite difference gradients for up to 50
# randomly selected reflections.
selection = random.sample(xrange(len(hkls)), min(len(hkls), 50))
uc = mycrystal.get_unit_cell()
exclusion_limit = max(uc.reciprocal_parameters()[0:3])

# prepare the prediction parameterisation with the scan-independent part of the
# current geometry
pred_param.prepare()

verbose = False
for iref in selection:

    hkl, angle, frame = hkls[iref], angles[iref], frames[iref]

    # re-predict this hkl based on the perturbed UB at its frame
    pred_param.compose(frame)

    UB = xlo_param.get_state() * xluc_param.get_state()

    ref_list = ref_predictor.predict(hkl, UB)

    if len(ref_list) == 0: continue

    if len(ref_list) == 1:
        angle = ref_list[0].rotation_angle
        s = matrix.col(ref_list[0].beam_vector)

    elif len(ref_list) == 2: # take the one with the closest angle
        phi_diff = (abs(ref_list[0].rotation_angle - angle),
                    abs(ref_list[1].rotation_angle - angle))
        if phi_diff[1] > phi_diff[0]:
            angle = ref_list[0].rotation_angle
            s = matrix.col(ref_list[0].beam_vector)
        else:
            angle = ref_list[1].rotation_angle
            s = matrix.col(ref_list[1].beam_vector)

    else: # cannot have more than two reflections, this should never execute
        raise RuntimeError("Predicted more than two angles for a single hkl")

    # get analytical gradients
    an_grads = pred_param.get_gradients(hkl, s, angle, panel_id, frame)

# NB, reflections that just touch the Ewald sphere have large
# derivatives of phi wrt some parameters (asymptotically approching
# infinity). Such reflections are likely to fail the comparison against
# FD gradients. Andrew Leslie points out that these reflections are not
# even integrated, because of susceptibility to errors in the Lorentz
# correction. To avoid failures of this test, we exclude reflections
# that are closer than the largest reciprocal lattice vector length from
# the plane of the rotation axis and beam direction.
    s0 = matrix.col(mybeam.get_s0())
    r = s - s0
    e_s0_plane_norm = matrix.col(
                    mygonio.get_rotation_axis()).cross(s0).normalize()
    r_dist_from_plane = abs(r.dot(e_s0_plane_norm))
    if r_dist_from_plane <= exclusion_limit:
        continue

    fd_grads = get_fd_gradients(pred_param, hkl, angle, frame, ref_predictor,
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
                print "Finite difference derivatives: %.6f, %.6f, %.6f" % tuple(fd_grad)
                finish_time = time()
                print "Time Taken: ",finish_time - start_time
                raise

finish_time = time()
print "Time Taken: ",finish_time - start_time

# if we got this far,
print "OK"

