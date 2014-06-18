"""We now have a set of DIALS models to replace my own classes

This script is an exploration of their interfaces in comparison with mine"""

from __future__ import division

from math import pi

# cctbx imports
from scitbx import matrix
from libtbx.test_utils import approx_equal

# import DIALS models
from dxtbx.model.experiment import Beam
from dxtbx.model.experiment import Panel
#from dxtbx.model.experiment import Detector
#from dxtbx.model.experiment import Goniometer

# but mostly use the factories instead for ease of use
from dxtbx.model.experiment import beam_factory, detector_factory, \
    goniometer_factory, scan_factory

# import my models
from rstbx.bpcx.detector_model.instrument_specifics import pilatus
from rstbx.bpcx import sensor
from dials.scratch.dgw.source_model import source
from dials.scratch.dgw.goniometer_model import goniometer

#####################################
# Compare instances of each version #
#####################################

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
# dials_beam has 'get_direction', 'get_s0', 'get_wavelength',
# 'get_unit_s0', 'set_unit_s0', 'set_direction', 'set_s0',
# 'set_wavelength'

# the analogue of my_beam.get_beam() is
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

print "OK"

##############################
# Compare goniometer methods #
##############################

# I only use the getter, get_axis. The equivalent in the dials models is
# get_rotation_axis.

assert my_gonio.get_axis() == matrix.col(dials_gonio.get_rotation_axis())

print "OK"

############################
# Compare detector methods #
############################

# I get a sensor out of my detector model like this
my_sensor = my_det.sensors()[0]

# the dials equivalent simply uses an array index
dials_panel = dials_det[0]

# direct use of the sensor constructor
lim = (0,50)
my_panel = sensor(origin, fast, slow, lim, lim)

# equivalent using the dials Panel
dials_panel = Panel("PAD", "Panel", fast, slow, origin,
            (lim[1]/200, lim[1]/200), (200, 200), (0, 2e20))

# The sensor model uses Python properties instead of getters.
assert my_panel.origin == dials_panel.get_origin()
assert my_panel.dir1 == dials_panel.get_fast_axis()
assert my_panel.dir2 == dials_panel.get_slow_axis()
assert my_panel.normal == dials_panel.get_normal()
assert my_panel.distance == dials_panel.get_distance()
assert my_panel.D == dials_panel.get_D_matrix()

# the frame limits are handled differently. The dials panel implicitly
# has a limit at the origin, whereas my sensor has an arbitrary letterbox
dp_lim = dials_panel.get_image_size_mm()
assert my_panel.lim1 == (0, dp_lim[0])
assert my_panel.lim2 == (0, dp_lim[1])

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
                      exposure_times = 0.1,
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

# actually we use the array index usually, so test that too
tst = myscan.get_array_index_from_angle(0.5, deg=False)
assert myscan.get_angle_from_array_index(tst, deg=False) == 0.5

# Note frame numbers start from 1.0, so the start of a frame is int n, its
# centre is n + 0.5 and the start of the next frame is n + 1.0
assert myscan.get_image_index_from_angle(0.0, deg=False) == 1.0

print "OK"
