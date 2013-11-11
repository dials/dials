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
from dials.algorithms.refinement.refinement_helpers import print_model_geometry

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
    """Provides get_UB(image_number) for scan-varying prediction"""

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
refs1 = ref_predictor.predict(indices)
refs2 = sv_predictor()

assert len(refs1) == len(refs2)
print "OK"

def sort_refs(reflections):

    """Sort reflections by Miller index and entering flag"""
    refs_sorted = sorted(reflections, key=lambda x: x.entering)
    refs_sorted = sorted(refs_sorted, key=lambda x: x.miller_index[2])
    refs_sorted = sorted(refs_sorted, key=lambda x: x.miller_index[1])
    refs_sorted = sorted(refs_sorted, key=lambda x: x.miller_index[0])

    return refs_sorted

refs1_sorted = sort_refs(refs1)
refs2_sorted = sort_refs(refs2)

for (r1, r2) in zip(refs1_sorted, refs2_sorted):
    assert r1.miller_index == r2.miller_index
    dphi = r1.rotation_angle - r2.rotation_angle
    assert abs(dphi) < 0.01 * im_width

print "OK"
