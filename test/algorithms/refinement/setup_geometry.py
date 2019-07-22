"""Setup experimental geometry for refinement test cases"""

# Python and cctbx imports
from __future__ import absolute_import, division, print_function

import random

from scitbx import matrix
from libtbx.phil import parse, command_line

# dxtbx experimental models
from dxtbx.model import BeamFactory
from dxtbx.model import GoniometerFactory
from dxtbx.model import DetectorFactory

# crystal model
from dxtbx.model import Crystal


# Local functions
def random_vector_close_to(vector, sd=0.5):
    return matrix.col(vector).rotate_around_origin(
        matrix.col((random.random(), random.random(), random.random())).normalize(),
        random.gauss(0, sd),
        deg=True,
    )


class Extract(object):
    """Parse and extract geometry model from PHIL"""

    def __init__(
        self, master_phil, local_overrides="", cmdline_args=None, verbose=False
    ):

        self._verbose = verbose

        arg_interpreter = command_line.argument_interpreter(master_phil=master_phil)

        user_phil = parse(local_overrides)
        cmdline_phils = []
        if cmdline_args:
            for arg in cmdline_args:
                cmdline_phils.append(arg_interpreter.process(arg))

        working_phil = master_phil.fetch(sources=[user_phil] + cmdline_phils)

        self._params = working_phil.extract().geometry.parameters

        self.set_seed()

        self.build_goniometer()

        self.build_crystal()

        self.build_beam()

        self.build_detector()

        # write changes back to the PHIL object
        temp = working_phil.extract()
        temp.geometry.parameters = self._params
        self.phil = master_phil.format(python_object=temp)

    def set_seed(self):

        if self._params.random_seed is not None:
            random.seed(self._params.random_seed)
            # set the flex random seed too
            from dials.array_family import flex

            flex.set_random_seed(self._params.random_seed)
            if self._verbose:
                msg = "Random seed set to %d while building models"
                print(msg % self._params.random_seed)

    def build_goniometer(self):

        self.goniometer = GoniometerFactory.known_axis(self._params.goniometer.axis)

    def build_beam(self):

        if self._params.beam.wavelength.random:
            wavelength = random.uniform(*self._params.beam.wavelength.range)
        else:
            wavelength = self._params.beam.wavelength.value

        assert self._params.beam.direction.method in [
            "inclination",
            "close_to",
            "exactly",
        ]

        if self._params.beam.direction.method == "inclination":

            if self._params.beam.direction.inclination.random:
                inclination = random.gauss(
                    0.0, self._params.beam.direction.inclination.angle
                )
            else:
                inclination = self._params.beam.direction.inclination.angle

            beam_dir = matrix.col((0, 0, 1)).rotate_around_origin(
                matrix.col((0, 1, 0)), inclination, deg=True
            )

        elif self._params.beam.direction.method == "close_to":

            temp = self._params.beam.direction.close_to.direction
            beam_dir = random_vector_close_to(
                temp, sd=self._params.beam.direction.close_to.sd
            )

        elif self._params.beam.direction.method == "exactly":

            beam_dir = matrix.col(self._params.beam.direction.exactly)

        self.beam = BeamFactory.make_beam(unit_s0=beam_dir, wavelength=wavelength)

    def build_detector(self):

        assert self._params.detector.directions.method in ["close_to", "exactly"]

        if self._params.detector.directions.method == "close_to":

            temp = self._params.detector.directions.close_to.dir1
            dir1 = random_vector_close_to(
                temp, sd=self._params.detector.directions.close_to.sd
            )

            n = random_vector_close_to(
                self._params.detector.directions.close_to.norm,
                sd=self._params.detector.directions.close_to.sd,
            )

        elif self._params.detector.directions.method == "exactly":

            temp = self._params.detector.directions.exactly.dir1
            dir1 = matrix.col(temp)

            n = matrix.col(self._params.detector.directions.exactly.norm)

        dir2 = n.cross(dir1).normalize()

        assert self._params.detector.centre.method in ["close_to", "exactly"]

        if self._params.detector.centre.method == "close_to":

            centre = random_vector_close_to(
                self._params.detector.centre.close_to.value,
                sd=self._params.detector.centre.close_to.sd,
            )

        elif self._params.detector.centre.method == "exactly":

            temp = self._params.detector.centre.exactly.value
            centre = matrix.col(temp)

        origin = centre - (
            0.5 * self._params.detector.npx_fast * self._params.detector.pix_size * dir1
            + 0.5
            * self._params.detector.npx_slow
            * self._params.detector.pix_size
            * dir2
        )
        self.detector = DetectorFactory.make_detector(
            "PAD",
            dir1,
            dir2,
            origin,
            (self._params.detector.pix_size, self._params.detector.pix_size),
            (self._params.detector.npx_fast, self._params.detector.npx_slow),
            (0, 1.0e6),
        )

    @staticmethod
    def _build_cell_vec(vec):

        if vec.length.random:
            length = random.uniform(*vec.length.range)
        else:
            length = vec.length.value

        assert vec.direction.method in ["close_to", "exactly"]

        if vec.direction.method == "close_to":

            x = random_vector_close_to(
                vec.direction.close_to.direction, sd=vec.direction.close_to.sd
            )

        elif vec.direction.method == "exactly":

            x = matrix.col(vec.direction.exactly.direction)

        return length * x

    def build_crystal(self):

        vecs = [
            self._build_cell_vec(axis)
            for axis in (
                self._params.crystal.a,
                self._params.crystal.b,
                self._params.crystal.c,
            )
        ]

        sg = self._params.crystal.space_group_symbol

        self.crystal = Crystal(*vecs, space_group_symbol=sg)
