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
from math import pi
from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.test_utils import approx_equal
from libtbx.phil import parse
from scitbx.array_family import flex

##### Import model builder

from setup_geometry import Extract

##### Imports for reflection prediction

from dials.algorithms.spot_prediction import IndexGenerator, ray_intersection
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction import ScansRayPredictor, \
  ExperimentsPredictor

#### Import model parameterisations

from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisation
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
s0_param = BeamParameterisation(mybeam, mygonio)
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

# Generate reflections
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                      space_group(space_group_symbols(1).hall()).type(),
                      resolution)
indices = index_generator.to_array()

# Predict rays within the sweep range
ray_predictor = ScansRayPredictor(experiments, sweep_range)
obs_refs = ray_predictor.predict(indices)

# Take only those rays that intersect the detector
intersects = ray_intersection(mydetector, obs_refs)
obs_refs = obs_refs.select(intersects)

# Make a reflection predictor and re-predict for all these reflections. The
# result is the same, but we gain also the flags and xyzcal.px columns
ref_predictor = ExperimentsPredictor(experiments)
obs_refs['id'] = flex.int(len(obs_refs), 0)
obs_refs = ref_predictor.predict(obs_refs)

# Set 'observed' centroids from the predicted ones
obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm']

# Invent some variances for the centroid positions of the simulated data
im_width = 0.1 * pi / 180.
px_size = mydetector[0].get_pixel_size()
var_x = flex.double(len(obs_refs), (px_size[0] / 2.)**2)
var_y = flex.double(len(obs_refs), (px_size[1] / 2.)**2)
var_phi = flex.double(len(obs_refs), (im_width / 2.)**2)
obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x, var_y, var_phi)

# use a ReflectionManager to exclude reflections too close to the spindle
from dials.algorithms.refinement.reflection_manager import ReflectionManager
refman = ReflectionManager(obs_refs, experiments, outlier_detector=None)

# Redefine the reflection predictor to use the type expected by the Target class
ref_predictor = ExperimentsPredictor(experiments)

# make a target to ensure reflections are predicted and refman is finalised
from dials.algorithms.refinement.target import \
  LeastSquaresPositionalResidualWithRmsdCutoff
target = LeastSquaresPositionalResidualWithRmsdCutoff(experiments,
    ref_predictor, refman, pred_param, restraints_parameterisation=None)

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
  assert approx_equal(x_grads, an_grads[i]["dX_dp"], eps=5.e-6)
  assert approx_equal(y_grads, an_grads[i]["dY_dp"], eps=5.e-6)
  assert approx_equal(phi_grads, an_grads[i]["dphi_dp"], eps=5.e-6)

# return to the initial state
pred_param.set_param_vals(p_vals)

finish_time = time()
print "Time Taken: ",finish_time - start_time

# if we got this far,
print "OK"
