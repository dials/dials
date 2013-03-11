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
from setup_geometry import extract

# Reflection prediction
from dials.scratch.dgw.prediction import angle_predictor, impact_predictor
from rstbx.diffraction import full_sphere_indices
from cctbx.sgtbx import space_group, space_group_symbols

# Imports for the target function
from dials.scratch.dgw.refinement.target import reflection_manager

args = sys.argv[1:]

master_phil = parse("""
include file geometry.params
include file minimiser.params
""", process_includes=True)

overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

models = extract(master_phil, local_overrides=overrides, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

#############################
# Generate some reflections #
#############################

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
indices = full_sphere_indices(
    unit_cell = mycrystal.get_unit_cell(),
    resolution_limit = resolution,
    space_group = space_group(space_group_symbols(1).hall()))

# Select those that are excited in a 30 degree sweep and get their angles
UB = mycrystal.get_U() * mycrystal.get_B()
ap = angle_predictor(mycrystal, mybeam, mygonio, resolution)
obs_indices, obs_angles = ap.observed_indices_and_angles_from_angle_range(
    phi_start_rad = 0.0, phi_end_rad = pi/6., indices = indices)

# Project positions on camera
ip = impact_predictor(mydetector, mygonio, mybeam, mycrystal)
hkls, d1s, d2s, angles, svecs = ip.predict(obs_indices.as_vec3_double(),
                                       obs_angles)

# Invent some uncertainties
im_width = 0.1 * pi / 180.
sigd1s = [mydetector.px_size_fast() / 2.] * len(hkls)
sigd2s = [mydetector.px_size_slow() / 2.] * len(hkls)
sigangles = [im_width / 2.] * len(hkls)

# Build list of observations in the reflection manager. This classifies each
# reflection as passing into or out of the Ewald sphere
refman = reflection_manager(hkls, svecs,
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
for hkl, v in refman._H.items():

    for i, e in enumerate(v.exiting):

        r = v.Sc[i] - s0
        r_orig = r.rotate(mygonio.get_axis(), -1., deg=True)

        # is it inside the Ewald sphere (i.e. exiting)?
        test = (s0 + r_orig).length() < s0.length()
        assert(e == test)

print "OK"
