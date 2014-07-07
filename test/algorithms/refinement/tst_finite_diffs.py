#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Test analytical calculation of gradients of the target function versus finite
difference calculations"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
import random
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Experimental model builder
from setup_geometry import Extract

# We will set up a mock scan and a mock experiment list
from dxtbx.model.scan import scan_factory
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ScansRayPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation

# Imports for the target function
from dials.algorithms.refinement.target import \
  LeastSquaresPositionalResidualWithRmsdCutoff
from dials.algorithms.refinement.reflection_manager import ReflectionManager

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
    include scope dials.test.algorithms.refinement.geometry_phil
    """, process_includes=True)

models = Extract(master_phil, overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

# Build a mock scan for a 180 degree sweep of 0.1 degree images
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,1800),
                      exposure_times = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(1800),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
assert sweep_range == (0., pi)
assert approx_equal(im_width, 0.1 * pi / 180.)

experiments = ExperimentList()
experiments.append(Experiment(
      beam=mybeam, detector=mydetector, goniometer=mygonio,
      scan=myscan, crystal=mycrystal, imageset=None))

###########################
# Parameterise the models #
###########################

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, mygonio)
xlo_param = CrystalOrientationParameterisation(mycrystal)
xluc_param = CrystalUnitCellParameterisation(mycrystal)

########################################################################
# Link model parameterisations together into a parameterisation of the #
# prediction equation                                                  #
########################################################################

pred_param = XYPhiPredictionParameterisation(experiments,
    [det_param], [s0_param], [xlo_param], [xluc_param])

################################
# Apply known parameter shifts #
################################

# shift detector by 0.2 mm each translation and 2 mrad each rotation
det_p_vals = det_param.get_param_vals()
p_vals = [a + b for a, b in zip(det_p_vals,
                                 [2.0, 2.0, 2.0, 2.0, 2.0, 2.0])]
det_param.set_param_vals(p_vals)

# shift beam by 2 mrad in one axis
s0_p_vals = s0_param.get_param_vals()
p_vals = list(s0_p_vals)
p_vals[1] += 2.0
s0_param.set_param_vals(p_vals)

# rotate crystal a bit (=2 mrad each rotation)
xlo_p_vals = xlo_param.get_param_vals()
p_vals = [a + b for a, b in zip(xlo_p_vals, [2.0, 2.0, 2.0])]
xlo_param.set_param_vals(p_vals)

#############################
# Generate some reflections #
#############################

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

# build a reflection predictor
ref_predictor = ScansRayPredictor(experiments, sweep_range)

obs_refs = ref_predictor.predict(indices)

# Invent some variances for the centroid positions of the simulated data
im_width = 0.1 * pi / 180.
px_size = mydetector[0].get_pixel_size()
var_x = (px_size[0] / 2.)**2
var_y = (px_size[1] / 2.)**2
var_phi = (im_width / 2.)**2

for ref in obs_refs:

  # calc and set the impact position, assuming all reflections
  # intersect panel 0.
  impacts = mydetector[0].get_ray_intersection(ref.beam_vector)
  ref.image_coord_mm = impacts

  # set the 'observed' centroids
  ref.centroid_position = ref.image_coord_mm + (ref.rotation_angle, )

  # set the centroid variance
  ref.centroid_variance = (var_x, var_y ,var_phi)

  # set the frame number, calculated from rotation angle
  ref.frame_number = myscan.get_image_index_from_angle(
      ref.rotation_angle, deg=False)

###############################
# Undo known parameter shifts #
###############################

s0_param.set_param_vals(s0_p_vals)
det_param.set_param_vals(det_p_vals)
xlo_param.set_param_vals(xlo_p_vals)

#####################################
# Select reflections for refinement #
#####################################

refman = ReflectionManager(obs_refs.to_table(centroid_is_mm=True), experiments)

##############################
# Set up the target function #
##############################

# Redefine the reflection predictor to use the type expected by the Target class
from dials.algorithms.refinement.prediction import ExperimentsPredictor
ref_predictor = ExperimentsPredictor(experiments)

mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(
    experiments, ref_predictor, refman, pred_param)

# get the functional and gradients
mytarget.predict()
L, dL_dp = mytarget.compute_functional_and_gradients()
try:
  curvs = mytarget.curvatures()
except AttributeError:
  curvs = None

####################################
# Do FD calculation for comparison #
####################################

# function for calculating finite difference gradients of the target function
def get_fd_gradients(target, pred_param, deltas):
  """Calculate centered finite difference gradients for each of the
  parameters of the target function.

  "deltas" must be a sequence of the same length as the parameter list, and
  contains the step size for the difference calculations for each parameter.
  """

  p_vals = pred_param.get_param_vals()
  assert len(deltas) == len(p_vals)
  fd_grad = []
  fd_curvs = []

  try:
    do_curvs = callable(target.curvatures)
  except AttributeError:
    do_curvs = False

  for i in range(len(deltas)):
    val = p_vals[i]

    sys.stdout.flush()
    p_vals[i] -= deltas[i] / 2.
    pred_param.set_param_vals(p_vals)
    target.predict()

    rev_state = target.compute_functional_and_gradients()

    p_vals[i] += deltas[i]
    pred_param.set_param_vals(p_vals)

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

  # return to the initial state
  pred_param.set_param_vals(p_vals)

  return fd_grad, fd_curvs

# test normalised differences between FD and analytical calculations
fdgrads = get_fd_gradients(mytarget, pred_param, [1.e-7] * len(pred_param))
diffs = [a - b for a, b in zip(dL_dp, fdgrads[0])]
norm_diffs = tuple([a / b for a, b in zip(diffs, fdgrads[0])])
for e in norm_diffs: assert abs(e) < 0.001 # check differences less than 0.1%
print "OK"

# test normalised differences between FD curvatures and analytical least
# squares approximation. We don't expect this to be especially close
if curvs:
  diffs = [a - b for a, b in zip(curvs, fdgrads[1])]
  norm_diffs = tuple([a / b for a, b in zip(diffs, fdgrads[1])])
  for e in norm_diffs: assert abs(e) < 0.1 # check differences less than 10%
print "OK"
