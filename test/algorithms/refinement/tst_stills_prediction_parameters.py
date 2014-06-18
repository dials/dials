#!/usr/bin/env python

#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
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

##### Import model builder

from setup_geometry import Extract

##### Imports for reflection prediction

from dials.algorithms.spot_prediction import IndexGenerator
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction import ScansRayPredictor

#### Import model parameterisations

from dials.algorithms.refinement.parameterisation.prediction_parameters_stills \
  import StillsPredictionParameterisation

from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, \
    CrystalUnitCellParameterisation

from time import time
start_time = time()

#### Create models

args = sys.argv[1:]

# build models, with a larger crystal than default in order to get enough
# reflections on the 'still' image
overrides = """
geometry.parameters.crystal.a.length.range=40 50;
geometry.parameters.crystal.b.length.range=40 50;
geometry.parameters.crystal.c.length.range=40 50;
geometry.parameters.random_seed = 42"""

master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    """, process_includes=True)

models = Extract(master_phil, overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

# Build a mock scan for a 72 degree sweep

from dxtbx.model.scan import scan_factory
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,1),
                      exposure_times = 0.1,
                      oscillation = (0, 3.0),
                      epochs = range(1),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)

#### Create parameterisations of these models

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, mygonio)
xlo_param = CrystalOrientationParameterisation(mycrystal)
xluc_param = CrystalUnitCellParameterisation(mycrystal)

# Create a scans ExperimentList, only for generating reflections
experiments = ExperimentList()
experiments.append(Experiment(
      beam=mybeam, detector=mydetector, goniometer=mygonio, scan=myscan,
      crystal=mycrystal, imageset=None))

# Create a stills ExperimentList
stills_experiments = ExperimentList()
stills_experiments.append(Experiment(
      beam=mybeam, detector=mydetector, crystal=mycrystal, imageset=None))

#### Unit tests

# Create a ScansRayPredictor
ref_predictor = ScansRayPredictor(experiments, sweep_range)

# Generate reflections
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                      space_group(space_group_symbols(1).hall()).type(),
                      resolution)
indices = index_generator.to_array()
ref_list = ref_predictor.predict(indices)

# make a reflection_table
reflections = ref_list.to_table(centroid_is_mm=True)

# Build a prediction parameterisation
pred_param = StillsPredictionParameterisation(stills_experiments,
               detector_parameterisations = [det_param],
               beam_parameterisations = [s0_param],
               xl_orientation_parameterisations = [xlo_param],
               xl_unit_cell_parameterisations = [xluc_param])

# use a ReflectionManager to exclude reflections too close to the spindle
from dials.algorithms.refinement.reflection_manager import StillsReflectionManager
refman = StillsReflectionManager(reflections, stills_experiments, iqr_multiplier=None)

# Redefine the reflection predictor to use the type expected by the Target class
from dials.algorithms.refinement.prediction import ExperimentsPredictor
ref_predictor = ExperimentsPredictor(stills_experiments)

# make a target to ensure reflections are predicted and refman is finalised
from dials.algorithms.refinement.target_stills import \
  LeastSquaresStillsResidualWithRmsdCutoff
target = LeastSquaresStillsResidualWithRmsdCutoff(stills_experiments,
    ref_predictor, refman, pred_param)

# keep only those reflections that pass inclusion criteria and have predictions
reflections = refman.get_matches()

# get analytical gradients
an_grads = pred_param.get_gradients(reflections)

# get finite difference gradients
p_vals = pred_param.get_param_vals()
deltas = [1.e-7] * len(p_vals)

from scitbx.array_family import flex
for i in range(len(deltas)):

  val = p_vals[i]

  p_vals[i] -= deltas[i] / 2.
  pred_param.set_param_vals(p_vals)

  ref_predictor.update()
  ref_predictor.predict(reflections)

  x, y, _ = reflections['xyzcal.mm'].deep_copy().parts()
  delpsi = reflections['delpsical.rad'].deep_copy()
  rev_state = flex.vec3_double(x, y, delpsi)

  p_vals[i] += deltas[i]
  pred_param.set_param_vals(p_vals)

  ref_predictor.update()
  ref_predictor.predict(reflections)

  x, y, _ = reflections['xyzcal.mm'].deep_copy().parts()
  delpsi = reflections['delpsical.rad'].deep_copy()
  fwd_state = flex.vec3_double(x, y, delpsi)
  p_vals[i] = val

  fd = (fwd_state - rev_state)
  x_grads, y_grads, delpsi_grads = fd.parts()
  x_grads /= deltas[i]
  y_grads /= deltas[i]
  delpsi_grads /= deltas[i]

  # compare FD with analytical calculations
  from scitbx.math import five_number_summary
  print "\n\nParameter", i
  grads = (x_grads, y_grads, delpsi_grads)

  for name, idx in zip(["X", "Y", "DeltaPsi"], (0, 1, 2)):
    print name
    rel_error = []
    abs_error = []
    for j, (a, b) in enumerate(zip(grads[idx], an_grads[idx][i])):
      abs_error.append(a - b)
      if abs(a + b) < 1.e-8:
        continue
      else:
        rel_error.append(2.*(a - b)/(a + b))
        #rel_error.append((a - b)/max(abs(a), abs(b)))
    if rel_error:
      print "summary of relative errors: %9.6f %9.6f %9.6f %9.6f %9.6f" % five_number_summary(rel_error)
    print "summary of absolute errors: %9.6f %9.6f %9.6f %9.6f %9.6f" % five_number_summary(abs_error)

  #assert approx_equal(x_grads, an_grads[0][i], eps=5.e-6)
  #assert approx_equal(y_grads, an_grads[1][i], eps=5.e-6)
  #assert approx_equal(delpsi_grads, an_grads[2][i], eps=5.e-6)

# return to the initial state
pred_param.set_param_vals(p_vals)

finish_time = time()
print "Time Taken: ",finish_time - start_time

# if we got this far,
print "OK"
