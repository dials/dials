#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Trivial check for whether classification of reflections as exiting or
entering the Ewald sphere is done the right way round"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from libtbx.phil import parse
from scitbx import matrix

# Get class to build experimental models
from setup_geometry import Extract

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.scratch.dgw.prediction import ReflectionPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# Imports for the target function
from dials.scratch.dgw.refinement.target import ReflectionManager

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

# Select those that are excited in a 30 degree sweep and get angles
UB = mycrystal.get_U() * mycrystal.get_B()
sweep_range = (0., pi/6.)
ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)

obs_refs = ref_predictor.predict(indices)

# Pull out reflection data as lists
temp = [(ref.miller_index, ref.rotation_angle,
         matrix.col(ref.beam_vector)) for ref in obs_refs]
hkls, angles, svecs = zip(*temp)

# Project positions on camera
# currently assume all reflections intersect panel 0
impacts = [mydetector[0].get_ray_intersection(
                        ref.beam_vector) for ref in obs_refs]
d1s, d2s = zip(*impacts)

# Invent some uncertainties
im_width = 0.1 * pi / 180.
px_size = mydetector.get_pixel_size()
sigd1s = [px_size[0] / 2.] * len(hkls)
sigd2s = [px_size[1] / 2.] * len(hkls)
sigangles = [im_width / 2.] * len(hkls)

# Build list of observations in the reflection manager. This classifies each
# reflection as passing into or out of the Ewald sphere
refman = ReflectionManager(hkls, svecs,
                        d1s, sigd1s,
                        d2s, sigd2s,
                        angles, sigangles,
                        mybeam, mygonio)

# Also update the reflection manager with the observation data as perfect
# predictions. We do this as the scattering vectors are only stored for
# the predicted reflections
refman.update_predictions(hkls, svecs, d1s, d2s, angles)

# for each reflection, reconstitute its relp vector and rotate it back by 1
# degree, so that it matches the originally generated reflection. Now test
# whether this vector lies inside or outside the Ewald sphere. If outside then
# the reflection is entering. If inside then the reflection is exiting.

s0 = matrix.col(mybeam.get_s0())
spindle = matrix.col(mygonio.get_rotation_axis())
for hkl, v in refman._obs_pred_pairs.items():

    for i, e in enumerate(v.exiting):

        r = v.Sc[i] - s0
        r_orig = r.rotate(spindle, -1., deg=True)

        # is it inside the Ewald sphere (i.e. exiting)?
        test = (s0 + r_orig).length() < s0.length()
        assert(e == test)

print "OK"
