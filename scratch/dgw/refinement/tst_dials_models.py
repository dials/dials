"""We now have a set of DIALS models to replace my own classes

This script is an exploration of their interfaces in comparison with mine"""

# Import DIALS models
from dials.model.experiment import Beam
#from dials.model.experiment import Panel, Detector
from dials.model.experiment import Goniometer

# use the factories instead for ease of use
from dials.model.experiment import beam_factory, detector_factory, \
    goniometer_factory

# Import my models
from rstbx.bpcx.detector_model.instrument_specifics import pilatus
from rstbx.bpcx import sensor
from dials.scratch.dgw.source_model import source
from dials.scratch.dgw.crystal_model import crystal
from dials.scratch.dgw.goniometer_model import goniometer

from scitbx import matrix

# beam
beam_direction = matrix.col((0., 0., 1.))
wavelength = 1.2

my_beam = source(beam_direction, wavelength)

bf = beam_factory()
dials_beam = bf.make_beam(beam_direction, wavelength)

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

#goniometer
spindle_axis = matrix.col((1., 0., 0.))

my_gonio = goniometer(spindle_axis)

gf = goniometer_factory()
dials_gonio = gf.make_goniometer(spindle_axis,
    fixed_rotation = (1, 0, 0, 0, 1, 0, 0, 0, 1))
