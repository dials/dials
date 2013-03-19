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
from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.test_utils import approx_equal
from libtbx.phil import parse

##### Import model builder

from setup_geometry import Extract

##### Imports for reflection prediction

from dials.algorithms.spot_prediction import IndexGenerator
from dials.scratch.dgw.prediction import ReflectionPredictor

#### Import model parameterisations

from dials.scratch.dgw.refinement.prediction_parameters import \
    DetectorSpacePredictionParameterisation
from dials.scratch.dgw.refinement.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.scratch.dgw.refinement.source_parameters import \
    BeamParameterisationOrientation
from dials.scratch.dgw.refinement.crystal_parameters import \
    CrystalOrientationParameterisation, \
    CrystalUnitCellParameterisation

#### Local functions

def print_grads(grad_list):
    for i, grad in enumerate(grad_list):
        print ("Param %02d. Gradients: "
               "%.5f, %.5f, %.5f" % ((i,) + tuple(grad)))

# Functions required for finite difference calculations

def get_state(det, hkl, angle, reflection_predictor):
    '''reflection prediction for the current state of the models'''

    # update reflection_predictor with latest geometry
    reflection_predictor.update()

    # predict for this hkl
    refs = reflection_predictor.predict(hkl)

    # select which is nearest the observed angle
    deltas = [abs(ref.rotation_angle - angle) for ref in refs]
    new_ref = refs[deltas.index(min(deltas))]

    # predict impact position
    impact = det[0].get_ray_intersection(new_ref.beam_vector)

    return matrix.col((impact) + (new_ref.rotation_angle,))

def get_fd_gradients(pred_param, hkl, phi, reflection_predictor,
                     deltas):
    '''Calculate centered finite difference gradients for each of the
    parameters of the prediction parameterisation object, for reflection
    hkl at angle phi.

    "deltas" must be a sequence of the same length as the parameter list,
    and contains the step size for the difference calculations for each
    parameter.'''

    gon = pred_param._gonio
    src = pred_param._beam
    det = pred_param._detector
    xl = pred_param._crystal
    rp = reflection_predictor

    p_vals = pred_param.get_p()
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):

        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.
        pred_param.set_p(p_vals)
        rev_state = get_state(det, hkl, phi, rp)

        p_vals[i] += deltas[i]
        pred_param.set_p(p_vals)
        fwd_state = get_state(det, hkl, phi, rp)

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

master_phil = parse("""
    include file geometry.params
    """, process_includes=True)

models = Extract(master_phil, overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

#### Create parameterisations of these models

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam)
xlo_param = CrystalOrientationParameterisation(mycrystal)
xluc_param = CrystalUnitCellParameterisation(mycrystal)

#### Unit tests

# Build a prediction parameterisation with a single detector model
pred_param = DetectorSpacePredictionParameterisation(mydetector,
             mybeam, mycrystal, mygonio, [det_param])

# Check the accessors
assert len(pred_param) == 6
for (a, b) in zip(pred_param.get_p(), det_param.get_p()):
    assert a==b

pred_param.set_p([100., 1.0, 1.0, 0., 0., 0.])
for (a, b) in zip(pred_param.get_p(), det_param.get_p()):
    assert a==b

# Build a full global parameterisation
pred_param = DetectorSpacePredictionParameterisation(
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
sweep_range = (0., pi/5.)
ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio,
                                    sweep_range)
ref_list = ref_predictor.predict(indices)

# Pull out lists of required reflection data
temp = [(ref.miller_index, ref.rotation_angle,
         matrix.col(ref.beam_vector)) for ref in ref_list]
hkls, angles, s_vecs = zip(*temp)

# Project positions on camera. Currently assuming all reflections
# intersect panel 0
impacts = [mydetector[0].get_ray_intersection(
                        ref.beam_vector) for ref in ref_list]
d1s, d2s = zip(*impacts)

# Test get_state for the first reflection
tmp = get_state(mydetector, hkls[0], angles[0], ref_predictor)
for (a, b) in zip(tmp, (d1s[0], d2s[0], angles[0])):
    assert a == b

# Compare analytical and finite difference gradients for up to 50
# randomly selected reflections.
selection = random.sample(xrange(len(hkls)), min(len(hkls), 50))
uc = mycrystal.get_unit_cell()
exclusion_limit = max(uc.reciprocal_parameters()[0:3])

verbose = False
for iref in selection:
    hkl, s, angle = hkls[iref], s_vecs[iref], angles[iref]

    # get analytical gradients
    an_grads = pred_param.get_gradients(hkl, s, angle)

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

    fd_grads = get_fd_gradients(pred_param, hkl, angle, ref_predictor,
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

finish_time = time()
print "Time Taken: ",finish_time - start_time

# if we got this far,
print "OK"

