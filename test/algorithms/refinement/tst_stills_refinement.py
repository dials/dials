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
A simple test of stills refinement using fake data (that is not really still).
We only attempt to refine the crystal. The beam and detector are known.

"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal

# Get modules to build models and minimiser using PHIL
import setup_geometry
import setup_minimiser

# We will set up a mock scan and a mock experiment list
from dxtbx.model.scan import scan_factory
from dials.model.experiment.experiment_list import ExperimentList, Experiment

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
from dials.algorithms.refinement.prediction import ScansRayPredictor
from dials.algorithms.spot_prediction import ray_intersection
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation

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
    include scope dials.test.algorithms.refinement.minimiser_phil
    """, process_includes=True)

# build models, with a larger crystal than default in order to get enough
# reflections on the 'still' image
param = """
geometry.parameters.crystal.a.length.range=40 50;
geometry.parameters.crystal.b.length.range=40 50;
geometry.parameters.crystal.c.length.range=40 50;
geometry.parameters.random_seed = 42"""
models = setup_geometry.Extract(master_phil, cmdline_args = args,
                        local_overrides=param)

crystal = models.crystal
mydetector = models.detector
mygonio = models.goniometer
mybeam = models.beam

# Build a mock scan for a 3 degree image, which will later be treated as a still
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,1),
                      exposure_times = 0.1,
                      oscillation = (0, 3.0),
                      epochs = range(1),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
assert approx_equal(im_width, 3.0 * pi / 180.)

# Build experiment lists
stills_experiments = ExperimentList()
stills_experiments.append(Experiment(
      beam=mybeam, detector=mydetector, crystal=crystal, imageset=None))
scans_experiments = ExperimentList()
scans_experiments.append(Experiment(
      beam=mybeam, detector=mydetector, crystal=crystal, goniometer = mygonio,
      scan=myscan, imageset=None))

##########################################################
# Parameterise the models (only for perturbing geometry) #
##########################################################

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisationOrientation(mybeam, goniometer=None)
xlo_param = CrystalOrientationParameterisation(crystal)
xluc_param = CrystalUnitCellParameterisation(crystal)

# Fix beam to the X-Z plane (imgCIF geometry)
s0_param.set_fixed([True, False])

################################
# Apply known parameter shifts #
################################

# rotate crystal a bit (=5 mrad each rotation)
xlo_p_vals = []
p_vals = xlo_param.get_param_vals()
xlo_p_vals.append(p_vals)
new_p_vals = [a + b for a, b in zip(p_vals, [5., 5., 5.])]
xlo_param.set_param_vals(new_p_vals)

# change unit cell a bit (=1.0 Angstrom length upsets, 0.5 degree of
# gamma angle)
xluc_p_vals = []
p_vals = xluc_param.get_param_vals()
xluc_p_vals.append(p_vals)
cell_params = crystal.get_unit_cell().parameters()
cell_params = [a + b for a, b in zip(cell_params, [1.0, 1.0, 1.0, 0.0,
                                                   0.0, 0.5])]
new_uc = unit_cell(cell_params)
newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
S = symmetrize_reduce_enlarge(crystal.get_space_group())
S.set_orientation(orientation=newB)
X = tuple([e * 1.e5 for e in S.forward_independent_parameters()])
xluc_param.set_param_vals(X)

#############################
# Generate some reflections #
#############################

#print "Reflections will be generated with the following geometry:"
#print mybeam
#print mydetector
#print crystal

# All indices in a 2.0 Angstrom sphere for crystal1
resolution = 2.0
index_generator = IndexGenerator(crystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

# Build a reflection predictor
ref_predictor = ScansRayPredictor(scans_experiments, sweep_range)

obs_refs = ref_predictor.predict(indices, experiment_id=0)

#print "Total number of reflections excited", len(obs_refs)

# Invent some variances for the centroid positions of the simulated data
im_width = 0.1 * pi / 180.
px_size = mydetector[0].get_pixel_size()
var_x = (px_size[0] / 2.)**2
var_y = (px_size[1] / 2.)**2
var_phi = (im_width / 2.)**2

obs_refs = ray_intersection(scans_experiments[0].detector, obs_refs)
for ref in obs_refs:

  # set the 'observed' centroids
  ref.centroid_position = ref.image_coord_mm + (ref.rotation_angle, )

  # set the centroid variance
  ref.centroid_variance = (var_x, var_y, var_phi)

  # set the frame number, calculated from rotation angle
  ref.frame_number = myscan.get_image_index_from_angle(
      ref.rotation_angle, deg=False)

  # ensure the crystal number is set to zero (should be by default)
  ref.crystal = 0

#print "Total number of observations made", len(obs_refs)

###############################
# Undo known parameter shifts #
###############################

xlo_param.set_param_vals(xlo_p_vals[0])
xluc_param.set_param_vals(xluc_p_vals[0])

#print "Refinement will start from the following geometry:"
#print mybeam
#print mydetector
#print crystal

# make a refiner
from dials.framework.registry import Registry
sysconfig = Registry().config()
params = sysconfig.params()

# overrides to fix beam and detector
params.refinement.parameterisation.beam.fix="all"
params.refinement.parameterisation.detector.fix="all"

# Change this to get a plot
do_plot = False
if do_plot: params.refinement.refinery.track_parameter_correlation=True

verbosity = 0
from dials.algorithms.refinement.refiner import RefinerFactory
refiner = RefinerFactory.from_parameters_data_experiments(params,
  obs_refs.to_table(centroid_is_mm=True), stills_experiments, verbosity=verbosity)

# run refinement
history = refiner.run()

if do_plot:
  plt = refiner.parameter_correlation_plot(len(history.parameter_correlation)-1)
  plt.show()

#print "Refinement has completed with the following geometry:"
#expts = refiner.get_experiments()
#for beam in expts.beams(): print beam
#for detector in expts.detectors(): print detector
#for crystal in  expts.crystals(): print crystal
