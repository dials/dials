#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""
Regression test for refinement of beam, detector and crystal orientation
parameters using generated reflection positions from ideal geometry.

"""

# Python and cctbx imports
from __future__ import division
from math import pi
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get modules to build models and minimiser using PHIL
import setup_geometry
import setup_minimiser

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
    DetectorSpacePredictionParameterisation

# Imports for the target function
from dials.algorithms.refinement.target import \
    LeastSquaresPositionalResidualWithRmsdCutoff, ReflectionManager

#############################
# Setup experimental models #
#############################

override = """geometry.parameters
{
  beam.wavelength.random=False
  beam.wavelength.value=1.0
  beam.direction.inclination.random=False
  crystal.a.length.random=False
  crystal.a.length.value=12.0
  crystal.a.direction.method=exactly
  crystal.a.direction.exactly.direction=1.0 0.002 -0.004
  crystal.b.length.random=False
  crystal.b.length.value=14.0
  crystal.b.direction.method=exactly
  crystal.b.direction.exactly.direction=-0.002 1.0 0.002
  crystal.c.length.random=False
  crystal.c.length.value=13.0
  crystal.c.direction.method=exactly
  crystal.c.direction.exactly.direction=0.002 -0.004 1.0
  detector.directions.method=exactly
  detector.directions.exactly.dir1=0.99 0.002 -0.004
  detector.directions.exactly.norm=0.002 -0.001 0.99
  detector.centre.method=exactly
  detector.centre.exactly.value=1.0 -0.5 199.0
}"""

master_phil = parse("""
include scope dials.test.algorithms.refinement.geometry_phil
include scope dials.test.algorithms.refinement.minimiser_phil
""", process_includes=True)

models = setup_geometry.Extract(master_phil, local_overrides=override,
                                verbose=False)

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

########################################################################
# Link model parameterisations together into a parameterisation of the #
# prediction equation                                                  #
########################################################################

pred_param = DetectorSpacePredictionParameterisation(
mydetector, mybeam, mycrystal, mygonio, [det_param], [s0_param],
[xlo_param], [xluc_param])

################################
# Apply known parameter shifts #
################################

# shift detector by 1.0 mm each translation and 4 mrad each rotation
det_p_vals = det_param.get_param_vals()
p_vals = [a + b for a, b in zip(det_p_vals,
                                [1.0, 1.0, 1.0, 4., 4., 4.])]
det_param.set_param_vals(p_vals)

# shift beam by 4 mrad in free axis
s0_p_vals = s0_param.get_param_vals()
p_vals = list(s0_p_vals)

p_vals[0] += 4.
s0_param.set_param_vals(p_vals)

# rotate crystal a bit (=3 mrad each rotation)
xlo_p_vals = xlo_param.get_param_vals()
p_vals = [a + b for a, b in zip(xlo_p_vals, [3., 3., 3.])]
xlo_param.set_param_vals(p_vals)

# change unit cell a bit (=0.1 Angstrom length upsets, 0.1 degree of
# alpha and beta angles)
xluc_p_vals = xluc_param.get_param_vals()
cell_params = mycrystal.get_unit_cell().parameters()
cell_params = [a + b for a, b in zip(cell_params, [0.1, -0.1, 0.1, 0.1,
                                                   -0.1, 0.0])]
new_uc = unit_cell(cell_params)
newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
S = symmetrize_reduce_enlarge(mycrystal.get_space_group())
S.set_orientation(orientation=newB)
X = tuple([e * 1.e5 for e in S.forward_independent_parameters()])
xluc_param.set_param_vals(X)

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
myscan = sf.make_scan(image_range = (1,1800),
                      exposure_time = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(1800),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
assert sweep_range == (0., pi)
assert approx_equal(im_width, 0.1 * pi / 180.)

ref_predictor = ReflectionPredictor(mycrystal, mybeam, mygonio, sweep_range)

obs_refs = ref_predictor.predict(indices)

# Pull out reflection data as lists
temp = [(ref.miller_index, ref.entering, ref.rotation_angle,
         matrix.col(ref.beam_vector)) for ref in obs_refs]
hkls, entering_flags, angles, svecs = zip(*temp)

# convert angles to image number
frames = map(lambda x: myscan.get_image_index_from_angle(x, deg=False),
             angles)

# Project positions on camera
# currently assume all reflections intersect panel 0
impacts = [mydetector[0].get_ray_intersection(
                        ref.beam_vector) for ref in obs_refs]
d1s, d2s = zip(*impacts)

print "Total number of observations made", len(hkls)

# Invent some uncertainties
im_width = 0.1 * pi / 180.
px_size = mydetector.get_pixel_size()
sigd1s = [px_size[0] / 2.] * len(hkls)
sigd2s = [px_size[1] / 2.] * len(hkls)
sigangles = [im_width / 2.] * len(hkls)

###############################
# Undo known parameter shifts #
###############################

s0_param.set_param_vals(s0_p_vals)
det_param.set_param_vals(det_p_vals)
xlo_param.set_param_vals(xlo_p_vals)
xluc_param.set_param_vals(xluc_p_vals)

#####################################
# Select reflections for refinement #
#####################################

refman = ReflectionManager(hkls, entering_flags, frames, svecs,
                           d1s, sigd1s,
                           d2s, sigd2s,
                           angles, sigangles,
                           mybeam, mygonio, myscan)

##############################
# Set up the target function #
##############################

# The current 'achieved' criterion compares RMSD against 1/3 the pixel size and
# 1/3 the image width in radians. For the simulated data, these are just made up
mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(ref_predictor,
    mydetector, refman, pred_param, im_width)

######################################
# Set up the LSTBX refinement engine #
######################################

overrides="""minimiser.parameters.engine=GaussNewtonIterations
minimiser.parameters.verbosity=0
minimiser.parameters.logfile=None"""
refiner = setup_minimiser.Extract(master_phil,
                                  mytarget,
                                  pred_param,
                                  local_overrides = overrides).refiner

refiner.run()

assert mytarget.achieved()
assert refiner.get_num_steps() == 1
assert approx_equal(mytarget.rmsds(), (0.005082333635071939,
                                       0.004212919308923247,
                                       8.976084121839415e-05))
print "OK"

###############################
# Undo known parameter shifts #
###############################

s0_param.set_param_vals(s0_p_vals)
det_param.set_param_vals(det_p_vals)
xlo_param.set_param_vals(xlo_p_vals)
xluc_param.set_param_vals(xluc_p_vals)

######################################################
# Set up the LBFGS with curvatures refinement engine #
######################################################

overrides="""minimiser.parameters.engine=LBFGScurvs
minimiser.parameters.verbosity=0
minimiser.parameters.logfile=None"""
refiner = setup_minimiser.Extract(master_phil,
                                  mytarget,
                                  pred_param,
                                  local_overrides = overrides).refiner

refiner.run()

assert mytarget.achieved()
assert refiner.get_num_steps() == 9
assert approx_equal(mytarget.rmsds(), (0.0564093680961448,
                                       0.03432636086464995,
                                       0.0003569908632122291))
print "OK"
