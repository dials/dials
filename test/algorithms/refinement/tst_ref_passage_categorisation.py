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

# Get class to build experimental models
from setup_geometry import Extract

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ReflectionPredictor
from cctbx.sgtbx import space_group, space_group_symbols
from dials.algorithms.spot_prediction import ray_intersection

# Imports for the target function
from dials.algorithms.refinement.target import ReflectionManager

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

# Select those that are excited in a 30 degree sweep and get angles
UB = mycrystal.get_U() * mycrystal.get_B()
sweep_range = (0., pi/6.)
ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)

obs_refs = ref_predictor.predict(indices)

# Project positions on camera
# currently assume all reflections intersect panel 0
impacts = ray_intersection(mydetector, obs_refs, panel=0)

# Pull out reflection data as lists
temp = [(ref.miller_index, ref.entering, ref.rotation_angle,
         ref.panel_number,
         matrix.col(ref.beam_vector)) for ref in obs_refs]
hkls, entering_flags, angles, panels, svecs = zip(*temp)

# Build a mock scan for a 30 degree sweep
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,300),
                      exposure_time = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(300),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
print sweep_range
assert approx_equal(sweep_range, (0., pi / 6.))
assert approx_equal(im_width, 0.1 * pi / 180.)

# convert angles to image number
frames = map(lambda x: myscan.get_image_index_from_angle(x, deg=False),
             angles)

# Pull out impact positions as lists
temp = [ref.image_coord_mm for ref in impacts]
d1s, d2s = zip(*temp)

# Invent some uncertainties
im_width = 0.1 * pi / 180.
px_size = mydetector.get_pixel_size()
sigd1s = [px_size[0] / 2.] * len(hkls)
sigd2s = [px_size[1] / 2.] * len(hkls)
sigangles = [im_width / 2.] * len(hkls)

# Build list of observations in the reflection manager. This classifies each
# reflection as passing into or out of the Ewald sphere
refman = ReflectionManager(hkls, entering_flags, frames, svecs,
                           panels,
                           d1s, sigd1s,
                           d2s, sigd2s,
                           angles, sigangles,
                           mybeam, mygonio, myscan)

mypanel = mydetector[0]
s0 = matrix.col(mybeam.get_s0())
spindle = matrix.col(mygonio.get_rotation_axis())

# for each reflection, reconstitute its relp vector and rotate it back by 1
# degree, so that it matches the originally generated reflection. Now test
# whether this vector lies inside or outside the Ewald sphere. If outside then
# the reflection is entering. If inside then the reflection is exiting.
for h in refman.get_indices():
    for obs in refman.get_obs(h):

        # get the s vector of this reflection
        tmp = matrix.col(mypanel.get_lab_coord((obs.Xo, obs.Yo))).normalize()
        s = tmp / mybeam.get_wavelength()

        r = s - s0
        r_orig = r.rotate(spindle, -1., deg=True)

        # is it outside the Ewald sphere (i.e. entering)?
        test = (s0 + r_orig).length() > s0.length()
        assert(obs.entering == test)

print "OK"
