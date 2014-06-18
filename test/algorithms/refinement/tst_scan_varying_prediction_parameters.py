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
from scitbx import matrix
from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.test_utils import approx_equal
from libtbx.phil import parse

##### Import model builders

from setup_geometry import Extract
from dxtbx.model.scan import scan_factory

##### Imports for reflection prediction

from dials.algorithms.spot_prediction import IndexGenerator
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction import ScansRayPredictor

#### Import model parameterisations

from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters import \
    VaryingCrystalPredictionParameterisation, VaryingCrystalPredictionParameterisationFast
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters import \
    ScanVaryingCrystalOrientationParameterisation, \
    ScanVaryingCrystalUnitCellParameterisation

#### Import helper functions

from dials.algorithms.refinement.refinement_helpers import random_param_shift

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

# Create parameterisations of these models, with 5 samples for the
# scan-varying crystal parameterisations

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, mygonio)
xlo_param = ScanVaryingCrystalOrientationParameterisation(
        mycrystal, myscan.get_array_range(), 5)
xluc_param = ScanVaryingCrystalUnitCellParameterisation(
        mycrystal, myscan.get_array_range(), 5)

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

# Generate an ExperimentList
experiments = ExperimentList()
experiments.append(Experiment(
      beam=mybeam, detector=mydetector, goniometer=mygonio, scan=myscan,
      crystal=mycrystal, imageset=None))
sweep_range = myscan.get_oscillation_range(deg=False)

# Create a ScansRayPredictor
ref_predictor = ScansRayPredictor(experiments, sweep_range)

#### Unit tests

# Build a prediction equation parameterisation
#pred_param = VaryingCrystalPredictionParameterisation(experiments, [det_param],
#                                        [s0_param], [xlo_param], [xluc_param])
# Use the 'fast' version as this is the default and is expected to have larger
# errors in the analytical gradients
pred_param = VaryingCrystalPredictionParameterisationFast(experiments, [det_param],
                                        [s0_param], [xlo_param], [xluc_param])

# Generate some reflections
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                      space_group(space_group_symbols(1).hall()).type(),
                      resolution)
indices = index_generator.to_array()
ref_list = ref_predictor.predict(indices)

# make a reflection_table
reflections = ref_list.to_table(centroid_is_mm=True)
# take 100 random reflections for speed
from scitbx.array_family import flex
reflections = reflections.select(flex.random_selection(len(reflections), 100))

# use a ReflectionManager to exclude reflections too close to the spindle,
# plus set the frame numbers
from dials.algorithms.refinement.reflection_manager import ReflectionManager
refman = ReflectionManager(reflections, experiments, iqr_multiplier=None)

# Redefine the reflection predictor to use the type expected by the Target class
from dials.algorithms.refinement.prediction import ExperimentsPredictor
ref_predictor = ExperimentsPredictor(experiments)

# make a target to ensure reflections are predicted and refman is finalised
from dials.algorithms.refinement.target import \
  LeastSquaresPositionalResidualWithRmsdCutoff
target = LeastSquaresPositionalResidualWithRmsdCutoff(experiments,
    ref_predictor, refman, pred_param)

# keep only those reflections that pass inclusion criteria and have predictions
reflections = refman.get_matches()

# get analytical gradients
pred_param.compose(reflections)
an_grads = pred_param.get_gradients(reflections)

# get finite difference gradients
p_vals = pred_param.get_param_vals()
deltas = [1.e-7] * len(p_vals)

for i in range(len(deltas)):

  val = p_vals[i]

  p_vals[i] -= deltas[i] / 2.
  pred_param.set_param_vals(p_vals)
  pred_param.compose(reflections)

  ref_predictor.update()
  ref_predictor.predict(reflections)

  rev_state = reflections['xyzcal.mm'].deep_copy()

  p_vals[i] += deltas[i]
  pred_param.set_param_vals(p_vals)
  pred_param.compose(reflections)

  ref_predictor.update()
  ref_predictor.predict(reflections)

  fwd_state = reflections['xyzcal.mm'].deep_copy()
  p_vals[i] = val

  fd = (fwd_state - rev_state)
  x_grads, y_grads, phi_grads = fd.parts()
  x_grads /= deltas[i]
  y_grads /= deltas[i]
  phi_grads /= deltas[i]

  for n, (a,b) in enumerate(zip(x_grads, an_grads[0][i])):
    assert approx_equal(a, b, eps=5.e-6)
  for n, (a,b) in enumerate(zip(y_grads, an_grads[1][i])):
    assert approx_equal(a, b, eps=5.e-6)
  for n, (a,b) in enumerate(zip(phi_grads, an_grads[2][i])):
    assert approx_equal(a, b, eps=5.e-6)

  # compare with analytical calculation
  #assert approx_equal(x_grads, an_grads[0][i], eps=5.e-6)
  #assert approx_equal(y_grads, an_grads[1][i], eps=5.e-6)
  #assert approx_equal(phi_grads, an_grads[2][i], eps=5.e-6)

# return to the initial state
pred_param.set_param_vals(p_vals)
pred_param.compose(reflections)

finish_time = time()
print "Time Taken: ",finish_time - start_time

# if we got this far,
print "OK"
