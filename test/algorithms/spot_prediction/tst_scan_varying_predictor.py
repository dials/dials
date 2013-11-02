#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test prediction of reflections using the scan-varying reflection
predictor.
"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get modules to build models and minimiser using PHIL
from dials.test.algorithms.refinement import setup_geometry

# We will set up a mock scan
from dxtbx.model.scan import scan_factory

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ReflectionPredictor
from dials.algorithms.refinement.prediction import \
    ScanVaryingReflectionListGenerator
from cctbx.sgtbx import space_group, space_group_symbols

# Imports for the target function
from dials.algorithms.refinement.target import \
    LeastSquaresPositionalResidualWithRmsdCutoff, ReflectionManager

# Import helper functions
from dials.algorithms.refinement import print_model_geometry

#############################
# Setup experimental models #
#############################

args = sys.argv[1:]
master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    """, process_includes=True)

models = setup_geometry.Extract(master_phil, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

class DummyPredictionParameterisation(object):
    '''Provides get_UB(image_number) for scan-varying prediction'''

    def __init__(self, crystal):

        self._crystal = crystal

    def get_UB(self, image_number):

        UB = matrix.sqr(self._crystal.get_U()) * \
             matrix.sqr(self._crystal.get_B())
        return UB

pred_param = DummyPredictionParameterisation(mycrystal)

#############################
# Generate some reflections #
#############################

print "Reflections will be generated with the following geometry:"
print_model_geometry(mybeam, mydetector, mycrystal)
print

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

# Build a mock scan for a 180 degree sweep
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,180),
                      exposure_time = 0.1,
                      oscillation = (0, 1.0),
                      epochs = range(180),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
assert sweep_range == (0., pi)
assert approx_equal(im_width, 1.0 * pi / 180.)

ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)

dmin = mydetector.get_max_resolution(mybeam.get_s0())
sv_predictor = ScanVaryingReflectionListGenerator(pred_param, mybeam,
                                            mygonio, myscan, resolution)
obs_refs = ref_predictor.predict(indices)

sv_obs_refs = sv_predictor()

print "Total number of reflections excited", len(obs_refs)
print "Total number of reflections excited by scan-varying predictor", len(sv_obs_refs)

# Invent some variances for the centroid positions of the simulated data
#im_width = 0.1 * pi / 180.
#px_size = mydetector[0].get_pixel_size()
#var_x = (px_size[0] / 2.)**2
#var_y = (px_size[1] / 2.)**2
#var_phi = (im_width / 2.)**2

#obs_refs = ray_intersection(mydetector, obs_refs)
#for ref in obs_refs:
#
#    # set the centroid variance
#    ref.centroid_variance = (var_x, var_y ,var_phi)
#
#    # set the frame number, calculated from rotation angle
#    ref.frame_number = myscan.get_image_index_from_angle(
#        ref.rotation_angle, deg=False)

#print "Total number of observations made", len(obs_refs)
