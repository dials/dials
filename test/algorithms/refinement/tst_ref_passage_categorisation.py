#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Trivial check for whether classification of reflections as exiting or
entering the Ewald sphere is done the right way round"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from libtbx.phil import parse
from scitbx import matrix
from libtbx.test_utils import approx_equal

# Building experimental models
from setup_geometry import Extract
from dials.model.experiment.experiment_list import ExperimentList, Experiment

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ScansRayPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# We will set up a mock scan
from dxtbx.model.scan import scan_factory

args = sys.argv[1:]

master_phil = parse("""
include scope dials.test.algorithms.refinement.geometry_phil
include scope dials.test.algorithms.refinement.minimiser_phil
""", process_includes=True)

overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

models = Extract(master_phil, local_overrides=overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

#############################
# Generate some reflections #
#############################

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

# Build a mock scan for a 30 degree sweep
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,300),
                      exposure_times = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(300),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
assert approx_equal(sweep_range, (0., pi / 6.))
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
assert approx_equal(im_width, 0.1 * pi / 180.)

# Create an ExperimentList for ScansRayPredictor
experiments = ExperimentList()
experiments.append(Experiment(
        beam=mybeam, detector=mydetector, goniometer=mygonio,
        scan=myscan, crystal=mycrystal, imageset=None))

# Select those that are excited in a 30 degree sweep and get angles
UB = mycrystal.get_U() * mycrystal.get_B()
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

print "Total number of observations made", len(obs_refs)

mypanel = mydetector[0]
s0 = matrix.col(mybeam.get_s0())
spindle = matrix.col(mygonio.get_rotation_axis())

for ref in obs_refs:

  # get the s vector of this reflection
  s = matrix.col(ref.beam_vector)

  r = s - s0
  r_orig = r.rotate(spindle, -1., deg=True)

  # is it outside the Ewald sphere (i.e. entering)?
  test = (s0 + r_orig).length() > s0.length()
  assert(ref.entering == test)

print "OK"
