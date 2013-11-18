#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Using a known, simple geometry, produce data for a plot of the derivatives of
reflection coordinates (X, Y, phi) with respect to any parameter"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
import random
from scitbx import matrix

# Experimental models
from dials.model.experiment import beam_factory
from dials.model.experiment import goniometer_factory
from dials.model.experiment import detector_factory
from dials.model.experiment.crystal_model import Crystal

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ReflectionPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation

# Import helper functions
from dials.algorithms.refinement.refinement_helpers import print_model_geometry

# Local functions
def random_direction_close_to(vector, sd = 0.5):
  return vector.rotate_around_origin(matrix.col(
              (random.random(),
               random.random(),
               random.random())).normalize(),
               random.gauss(0, sd),  deg = True)

##############################
# Set and report random seed #
##############################

if len(sys.argv) == 1:
  seed = random.randint(0, sys.maxint)
elif len(sys.argv) == 2:
  seed = int(sys.argv[1])
else:
  print "Usage:", sys.argv[0], "[random_seed]"
  sys.exit(1)

random.seed(seed)
print "Random seed for this run:", seed

#############################
# Setup experimental models #
#############################

# Make a spindle along the X axis (imgCIF frame)
mygonio = goniometer_factory.known_axis((1, 0, 0))

# Put the beam along the Z axis
mybeam = beam_factory.make_beam(unit_s0=(0., 0., 1.), wavelength=1.0)

# Make a detector modelled on PILATUS 2M, S/N 24-0107 Diamond, 200 mm from the
# lab frame origin along the Z axis, with plane directions aligned with X and -Y axes
class dumb_sensor: pass
ds = dumb_sensor()
ds.dir1 = matrix.col((1., 0., 0.))
n = matrix.col((0., 0., 1.))
ds.dir2 = n.cross(ds.dir1).normalize()
ds.npx_fast = 1475
ds.npx_slow = 1679
ds.pix_size = 0.172
ds.centre = matrix.col((0., 0., 200.))
ds.origin = ds.centre - (0.5 * ds.npx_fast * ds.pix_size * ds.dir1 +
                         0.5 * ds.npx_slow * ds.pix_size * ds.dir2)
mydetector = detector_factory.make_detector("PAD",
                            ds.dir1, ds.dir2, ds.origin,
                            (ds.pix_size,ds.pix_size),
                            (ds.npx_fast, ds.npx_slow),
                            (0, 1.e6))

# Make a random P1 crystal
a = random.uniform(10,30) * random_direction_close_to(matrix.col((1, 0, 0)))
b = random.uniform(10,30) * random_direction_close_to(matrix.col((0, 1, 0)))
c = random.uniform(10,30) * random_direction_close_to(matrix.col((0, 0, 1)))
mycrystal = Crystal(a, b, c, space_group_symbol="P 1")

print "Reflections will be generated with the following geometry:"
print_model_geometry(mybeam, mydetector, mycrystal)

#############################
# Generate some reflections #
#############################

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

# Select those that are excited in a 90 degree sweep and get angles
UB = mycrystal.get_U() * mycrystal.get_B()
sweep_range = (0., pi/2.)
ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)

predictions = ref_predictor.predict(indices)

print "Total number of reflections excited", len(predictions)

# Project positions on camera
from dials.algorithms.spot_prediction import ray_intersection
obs_refs = ray_intersection(mydetector, predictions)

print "Total number of observations made", len(obs_refs)

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

pred_param = XYPhiPredictionParameterisation(
    mydetector, mybeam, mycrystal, mygonio, [det_param], [s0_param],
    [xlo_param], [xluc_param])

##########################################################################
# Get the gradients of X, Y, phi for each parameter, for each reflection #
##########################################################################

nparam = len(pred_param)

grads = [pred_param.get_gradients(ref.miller_index,
                                  ref.beam_vector,
                                  ref.rotation_angle,
                                  ref.panel_number) for ref in obs_refs]

f = open("plot_derivatives.dat", "w")
s1 = ["dX_dParam%02d"] * nparam
s2 = ["dY_dParam%02d"] * nparam
s3 = ["dPhi_dParam%02d"] * nparam
header_line = (("X, Y, Phi, " + ", ".join(s1 + s2 + s3) + "\n" ) % tuple(range(nparam) * 3))
f.write(header_line)
fmt = ", ".join(["%.5f"] * (3 + 3 * nparam)) + "\n"
#for x, y, phi, grad in zip(d1s, d2s, angles, grads):
for ref, grad in zip(obs_refs, grads):
  x = ref.image_coord_mm[0]
  y = ref.image_coord_mm[1]
  phi = ref.rotation_angle
  xgrads = [g[0] for g in grad]
  ygrads = [g[1] for g in grad]
  phigrads = [g[2] for g in grad]
  line = fmt % tuple([x] + [y] + [phi] + xgrads + ygrads + phigrads)
  f.write(line)
f.close()
