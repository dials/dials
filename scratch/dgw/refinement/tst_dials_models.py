"""We now have a set of DIALS models to replace my own classes

This script is an exploration of their interfaces in comparison with mine"""

from __future__ import division

from math import pi

# cctbx imports
from scitbx import matrix
from libtbx.test_utils import approx_equal

# import DIALS models
from dials.model.experiment import Beam
from dials.model.experiment import Panel
#from dials.model.experiment import Detector
#from dials.model.experiment import Goniometer

# but use the factories instead for ease of use
from dials.model.experiment import beam_factory, detector_factory, \
    goniometer_factory, scan_factory

# import my models
from rstbx.bpcx.detector_model.instrument_specifics import pilatus
from rstbx.bpcx import sensor
from dials.scratch.dgw.source_model import source
#from dials.scratch.dgw.crystal_model import Crystal
from dials.scratch.dgw.goniometer_model import goniometer

#############################
# Instances of each version #
#############################

# beam
beam_direction = matrix.col((0., 0., 1.))
sample_to_source = -1. * beam_direction
wavelength = 1.2

my_beam = source(beam_direction, wavelength)

# dials beam factory takes sample_to_source rather than beam_direction
# although we could use named arguments with unit_s0=beam_direction,
# wavelength=wavelength
bf = beam_factory()
dials_beam = bf.make_beam(sample_to_source, wavelength)

# detector
npx_fast = 1475
npx_slow = 1679
pix_size_f = pix_size_s = 0.172
fast = matrix.col((1., 0., 0.))
norm = matrix.col((0., 0., 1.))
slow = norm.cross(fast).normalize()
centre = matrix.col((0., 0., 200.))
origin = centre - (0.5 * npx_fast * pix_size_f * fast +
                   0.5 * npx_slow * pix_size_s * slow)

my_det = pilatus(origin, fast, slow, pix_size_f, pix_size_s, npx_fast, \
                        npx_slow)

df = detector_factory()
dials_det = df.make_detector("PAD", fast, slow, origin,
            (pix_size_f, pix_size_s), (npx_fast, npx_slow), (0, 2e20))

# goniometer
spindle_axis = matrix.col((1., 0., 0.))

my_gonio = goniometer(spindle_axis)

gf = goniometer_factory()
dials_gonio = gf.make_goniometer(spindle_axis,
    fixed_rotation = (1, 0, 0, 0, 1, 0, 0, 0, 1))

########################
# Compare beam methods #
########################

# my_beam has 'get_beam', 'get_s0', 'get_wavelength', 'set_s0' as its
# public interface.
# dials_beam has 'get_direction', 'get_s0', 'get_wavelength', 'get_unit_s0'
# 'set_unit_s0', 'set_direction', 'set_s0', 'set_wavelength'

# the analogue of my_beam.get_beam() is
print my_beam.get_beam()
print dials_beam.get_wavelength() * matrix.col(dials_beam.get_unit_s0())

assert my_beam.get_beam() == dials_beam.get_wavelength() * \
                             matrix.col(dials_beam.get_unit_s0())

# the dials model getters return a tuple, whereas my models return
# matrix.col
assert my_beam.get_s0() == matrix.col(dials_beam.get_s0())

# get_wavelength works the same way
assert my_beam.get_wavelength() == dials_beam.get_wavelength()

# in tests I directly use the constructor of source. It works too for
# the dials model
random_vec = matrix.col.random(3, 0.5, 1.5)
my_s0 = source(random_vec)
dials_s0 = Beam(random_vec)
assert approx_equal(my_s0.get_s0(), matrix.col(dials_s0.get_s0()))

# compare the setters. My set_s0 setter does not allow setting the
# wavelength, whereas the dials model does. To get the same behaviour as
# mine, use set_unit_s0 instead
random_vec = matrix.col.random(3, 0.5, 1.5)
my_s0.set_s0(random_vec)
dials_s0.set_unit_s0(random_vec)
assert approx_equal(my_s0.get_s0(), matrix.col(dials_s0.get_s0()))

# To change to DIALS models, I must modify these modules:
## prediction/tst
## prediction/predictors
## beam_parameters
## prediction_parameters
## tst_prediction_parameters
## setup_geometry
## refinement
## target
## plot_derivatives
## tst_finite_diffs
## tst_orientation_refinement
## tst_convergence_radius
## tst_convergence_radius_one_parameter
## tst_ref_passage_categorisation

# source_model is then retired, and probably should be removed

##############################
# Compare goniometer methods #
##############################

# I only use the getter, get_axis. The equivalent in the dials models is
# get_rotation_axis.

assert my_gonio.get_axis() == matrix.col(dials_gonio.get_rotation_axis())

# To change to DIALS models, I must modify these modules:
## prediction/tst
## prediction/predictors
## prediction_parameters
## tst_prediction_parameters
## setup_geometry
## target
## plot_derivatives

# tst_ref_passage_categorisation

############################
# Compare detector methods #
############################

# To change to DIALS models, I must modify these modules:

# prediction/predictors: change ImpactPredictor to use the prediction
# method built into the Panel object, that is get_ray_intersection(ray).
# Okay, in fact the whole prediction module should be replaced with
# DIALS code. There is already dials.algorithms.spot_prediction. I
# should look into this before continuing.

# Tests are written into prediction/tsts.py.
# It appears to behave as I want, so go ahead and use it. First convert
# tst_orientation refinement (fix the other tests another day).

# Now in my predictors.py I made a ReflectionPredictor class. This wraps
# DIALS RayPredictor, instantiating a new one when requested with
# updated experimental geometry. Use this class in Target to generate
# new predictions, and in tst_orientation_refinement to generate the
# 'observations'

# Next change detector_parameters then prediction_parameters to use
# DIALS classes, and test with tst_prediction_parameters.

# Changing tst_prediction_parameters requires changes to setup_geometry
# and conversion of the get_state function to use ReflectionPredictor
# rather than rstbx's reflection_prediction.

# Target also now needs to take a Detector, not an ImpactPredictor. So
# tst_orientation_refinement now will not work, until it also
# understands Detectors. Things to do:

# 1. fix setup_geometry to return a Detector
# 2. fix tst_prediction_parameters so that get_state uses
#   ReflectionPredictor from predictors and a Detector to calculate the
#   impacts. Also remove AnglePredictor_rstbx and replace with
#   ReflectionPredictor
# 3. change tst_orientation_refinement to use ReflectionPredictor too

# Orientation refinement now proceeds with DIALS classes.

# Next tidy up. I expect the following to be broken and require change:
# tst_finite_diffs FIXED
# tst_convergence_radius FIXED
# tst_convergence_radius_one_parameter FIXED
# tst_ref_passage_categorisation FIXED
# plot_derivatives

# FIXME the current detector parameterisation is initialised with a
# sensor, not a detector! This should be changed now, after I move to
# the DIALS models.

# I get a sensor out of my detector model like this
my_sensor = my_det.sensors()[0]

# the dials equivalent simply uses an array index
dials_panel = dials_det[0]

# direct use of the sensor constructor
lim = (0,50)
my_panel = sensor(origin, fast, slow, lim, lim)

# equivalent using the dials Panel
dials_panel = Panel("PAD", fast, slow, origin,
            (lim[1]/200, lim[1]/200), (200, 200), (0, 2e20))

# The sensor model uses Python properties instead of getters.
my_panel.origin == dials_panel.get_origin()
my_panel.dir1 == dials_panel.get_fast_axis()
my_panel.dir2 == dials_panel.get_slow_axis()
my_panel.normal == dials_panel.get_normal()
my_panel.distance == dials_panel.get_distance()
my_panel.D == dials_panel.get_D_matrix()

# the frame limits are handled differently. The dials panel implicitly
# has a limit at the origin, whereas my sensor has an arbitrary letterbox
dp_lim = dials_panel.get_image_size_mm()
my_panel.lim1 == (0, dp_lim[0])
my_panel.lim2 == (0, dp_lim[1])


# use of the sensor/Panel setter
my_panel.set_frame(origin, fast, slow)

dials_panel.set_frame(fast, slow, origin)

print "OK"

###################
# Test DIALS Scan #
###################

# Build a mock scan for a 180 degree sweep of 0.1 degree images
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

# check reciprocity of angle<-->frame functions
tst = myscan.get_image_index_from_angle(0.5, deg=False)
assert myscan.get_angle_from_image_index(tst, deg=False) == 0.5

# Note frame numbers start from 1.0, so the start of a frame is int n, its
# centre is n + 0.5 and the start of the next frame is n + 1.0
assert myscan.get_image_index_from_angle(0.0, deg=False) == 1.0

print "OK"
