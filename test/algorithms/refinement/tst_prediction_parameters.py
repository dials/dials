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
from dials.model.experiment.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction import ScansRayPredictor

#### Import model parameterisations

from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation
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

# Build a mock scan for a 72 degree sweep
sweep_range = (0., pi/5.)
from dxtbx.model.scan import scan_factory
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,720),
                      exposure_times = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(720),
                      deg = True)

#### Create parameterisations of these models

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, mygonio)
xlo_param = CrystalOrientationParameterisation(mycrystal)
xluc_param = CrystalUnitCellParameterisation(mycrystal)

# Create an ExperimentList
experiments = ExperimentList()
experiments.append(Experiment(
      beam=mybeam, detector=mydetector, goniometer=mygonio, scan=myscan,
      crystal=mycrystal, imageset=None))

#### Unit tests

# Build a prediction parameterisation
pred_param = XYPhiPredictionParameterisation(experiments,
               detector_parameterisations = [det_param],
               beam_parameterisations = [s0_param],
               xl_orientation_parameterisations = [xlo_param],
               xl_unit_cell_parameterisations = [xluc_param])

# Create a ScansRayPredictor
ref_predictor = ScansRayPredictor(experiments, sweep_range)

# Rebuild a full global parameterisation
pred_param = XYPhiPredictionParameterisation(experiments, [det_param],
                                        [s0_param], [xlo_param], [xluc_param])

# Generate reflections
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                      space_group(space_group_symbols(1).hall()).type(),
                      resolution)
indices = index_generator.to_array()
ref_list = ref_predictor.predict(indices)

# make a reflection_table
reflections = ref_list.to_table(centroid_is_mm=True)

# use a ReflectionManager to exclude reflections too close to the spindle
from dials.algorithms.refinement.reflection_manager import ReflectionManager
refman = ReflectionManager(reflections, experiments, iqr_multiplier=None)

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
an_grads = pred_param.get_gradients(reflections)

# get finite difference gradients
p_vals = pred_param.get_param_vals()
deltas = [1.e-7] * len(p_vals)

for i in range(len(deltas)):

  val = p_vals[i]

  p_vals[i] -= deltas[i] / 2.
  pred_param.set_param_vals(p_vals)

  ref_predictor.update()
  ref_predictor.predict(reflections)

  rev_state = reflections['xyzcal.mm'].deep_copy()
  #x_calc, y_calc, phi_calc = rev_state.parts()

  p_vals[i] += deltas[i]
  pred_param.set_param_vals(p_vals)

  ref_predictor.update()
  ref_predictor.predict(reflections)

  fwd_state = reflections['xyzcal.mm'].deep_copy()
  p_vals[i] = val

  fd = (fwd_state - rev_state)
  x_grads, y_grads, phi_grads = fd.parts()
  x_grads /= deltas[i]
  y_grads /= deltas[i]
  phi_grads /= deltas[i]

  # compare with analytical calculation
  assert approx_equal(x_grads, an_grads[0][i], eps=5.e-6)
  assert approx_equal(y_grads, an_grads[1][i], eps=5.e-6)
  assert approx_equal(phi_grads, an_grads[2][i], eps=5.e-6)

# return to the initial state
pred_param.set_param_vals(p_vals)

finish_time = time()
print "Time Taken: ",finish_time - start_time

# if we got this far,
print "OK"
