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
from scitbx.array_family import flex

# Building experimental models
from setup_geometry import Extract
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ScansRayPredictor, \
  ExperimentsPredictor
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
im_width = myscan.get_oscillation(deg=False)[1]
assert approx_equal(im_width, 0.1 * pi / 180.)

# Create an ExperimentList for ScansRayPredictor
experiments = ExperimentList()
experiments.append(Experiment(
        beam=mybeam, detector=mydetector, goniometer=mygonio,
        scan=myscan, crystal=mycrystal, imageset=None))

# Select those that are excited in a 30 degree sweep and get angles
UB = mycrystal.get_U() * mycrystal.get_B()
ray_predictor = ScansRayPredictor(experiments, sweep_range)
obs_refs = ray_predictor(indices)

# Set the experiment number
obs_refs['id'] = flex.int(len(obs_refs), 0)

# Calculate intersections
ref_predictor = ExperimentsPredictor(experiments)
obs_refs = ref_predictor(obs_refs)

print "Total number of observations made", len(obs_refs)

mypanel = mydetector[0]
s0 = matrix.col(mybeam.get_s0())
spindle = matrix.col(mygonio.get_rotation_axis())

for ref in obs_refs:

  # get the s1 vector of this reflection
  s1 = matrix.col(ref['s1'])

  r = s1 - s0
  r_orig = r.rotate(spindle, -1., deg=True)

  # is it outside the Ewald sphere (i.e. entering)?
  test = (s0 + r_orig).length() > s0.length()
  assert(ref['entering'] == test)

print "OK"
