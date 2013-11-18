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
Test refinement of beam, detector and crystal orientation parameters
using generated reflection positions from ideal geometry.

Control of the experimental model and choice of minimiser is done via
PHIL, which means we can do, for example:

cctbx.python tst_orientation_refinement.py \
"random_seed=3; engine=LBFGScurvs"

"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get modules to build models using PHIL
import setup_geometry

# We will set up a mock scan
from dxtbx.model.scan import scan_factory

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisationOrientation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints.parameter_reduction import \
    symmetrize_reduce_enlarge

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ReflectionPredictor
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation

# Imports for the target function
from dials.algorithms.refinement.target import \
    LeastSquaresPositionalResidualWithRmsdCutoff, ReflectionManager

# Import helper functions
from dials.algorithms.refinement.refinement_helpers \
    import print_model_geometry, refine

#############################
# Setup experimental models #
#############################

args = sys.argv[1:]
master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    include scope dials.test.algorithms.refinement.minimiser_phil
    """, process_includes=True)

models = setup_geometry.Extract(master_phil, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

###########################
# Parameterise the models #
###########################

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, mygonio)
xlo_param = CrystalOrientationParameterisation(mycrystal)
xluc_param = CrystalUnitCellParameterisation(mycrystal)

# Fix beam to the X-Z plane (imgCIF geometry)
s0_param.set_fixed([True, False])

# Fix crystal parameters
#xluc_param.set_fixed([True, True, True, True, True, True])

########################################################################
# Link model parameterisations together into a parameterisation of the #
# prediction equation                                                  #
########################################################################

pred_param = XYPhiPredictionParameterisation(
mydetector, mybeam, mycrystal, mygonio, [det_param], [s0_param],
[xlo_param], [xluc_param])

################################
# Apply known parameter shifts #
################################

# shift detector by 1.0 mm each translation and 2 mrad each rotation
det_p_vals = det_param.get_param_vals()
p_vals = [a + b for a, b in zip(det_p_vals,
                                [1.0, 1.0, 1.0, 2., 2., 2.])]
det_param.set_param_vals(p_vals)

# shift beam by 2 mrad in free axis
s0_p_vals = s0_param.get_param_vals()
p_vals = list(s0_p_vals)

p_vals[0] += 2.
s0_param.set_param_vals(p_vals)

# rotate crystal a bit (=2 mrad each rotation)
xlo_p_vals = xlo_param.get_param_vals()
p_vals = [a + b for a, b in zip(xlo_p_vals, [2., 2., 2.])]
xlo_param.set_param_vals(p_vals)

# change unit cell a bit (=0.1 Angstrom length upsets, 0.1 degree of
# gamma angle)
xluc_p_vals = xluc_param.get_param_vals()
cell_params = mycrystal.get_unit_cell().parameters()
cell_params = [a + b for a, b in zip(cell_params, [0.1, 0.1, 0.1, 0.0,
                                                   0.0, 0.1])]
new_uc = unit_cell(cell_params)
newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
S = symmetrize_reduce_enlarge(mycrystal.get_space_group())
S.set_orientation(orientation=newB)
X = tuple([e * 1.e5 for e in S.forward_independent_parameters()])
xluc_param.set_param_vals(X)

#############################
# Generate some reflections #
#############################

print "Reflections will be generated with the following geometry:"
print_model_geometry(mybeam, mydetector, mycrystal)
print "Target values of parameters are"
msg = "Parameters: " + "%.5f " * len(pred_param)
print msg % tuple(pred_param.get_param_vals())
print

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

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

ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)

excited_refs = ref_predictor.predict(indices)

print "Total number of reflections excited", len(excited_refs)

# Project positions on camera (assume panel 0) and set made up
# centroid variances and frame numbers
px_size = mydetector[0].get_pixel_size()

def set_impact(ref):
  """helper function to set centroid and fake variance and frame
  number in a reflection"""
  try:
    ref.image_coord_mm = mydetector[0].get_ray_intersection(ref.beam_vector)
    ref.centroid_position = ref.image_coord_mm + (ref.rotation_angle, )
    ref.centroid_variance = (px_size[0] / 2., px_size[1] / 2., im_width / 2.)
    ref.frame_number = myscan.get_image_index_from_angle(ref.rotation_angle, deg=False)
    return ref
  except RuntimeError:
    # this reflection doesn't intersect the Panel
    return False

obs_refs = [e for e in excited_refs if set_impact(e)]

print "Total number of observations made", len(obs_refs)

###############################
# Undo known parameter shifts #
###############################

s0_param.set_param_vals(s0_p_vals)
det_param.set_param_vals(det_p_vals)
xlo_param.set_param_vals(xlo_p_vals)
xluc_param.set_param_vals(xluc_p_vals)

print "Initial values of parameters are"
msg = "Parameters: " + "%.5f " * len(pred_param)
print msg % tuple(pred_param.get_param_vals())
print

###########################
# Use the refine function #
###########################

print "Prior to refinement the experimental model is:"
print_model_geometry(mybeam, mydetector, mycrystal)

refine(mybeam, mygonio, mycrystal, mydetector, myscan, obs_refs,
       verbosity=1, fix_cell=False, scan_varying=False)

print
print "Refinement has completed with the following geometry:"
print_model_geometry(mybeam, mydetector, mycrystal)
