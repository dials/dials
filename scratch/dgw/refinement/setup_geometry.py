#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Setup experimental geometry for refinement test cases"""

# Python and cctbx imports
from __future__ import division
import os, sys
import random
from scitbx import matrix
from libtbx.phil import parse, command_line

# Experimental models
from rstbx.bpcx.detector_model.instrument_specifics import pilatus
from dials.scratch.dgw.source_model import source
from dials.scratch.dgw.crystal_model import crystal
from dials.scratch.dgw.goniometer_model import goniometer

# Local functions
def random_direction_close_to(vector, sd = 0.5):
    return vector.rotate_around_origin(matrix.col(
                (random.random(),
                 random.random(),
                 random.random())).normalize(),
                 random.gauss(0, sd),  deg = True)

class extract(object):
    '''Parse and extract geometry model from PHIL'''

    def __init__(self, master_phil, local_overrides = "",
                 cmdline_args = None, verbose=True):

        self._verbose = verbose

        arg_interpreter = command_line.argument_interpreter(
            master_phil=master_phil)

        user_phil = parse(local_overrides)
        cmdline_phils = []
        if cmdline_args:
            for arg in cmdline_args:
                cmdline_phils.append(arg_interpreter.process(arg))

        working_phil = master_phil.fetch(
            sources=[user_phil] + cmdline_phils)

        self._params = working_phil.extract().geometry.parameters

        self.set_seed()

        self.build_goniometer()

        self.build_crystal()

        self.build_source()

        self.build_detector()

        # write changes back to the PHIL object
        temp = working_phil.extract()
        temp.geometry.parameters = self._params
        self.phil = master_phil.format(python_object = temp)

    def set_seed(self):

        if not self._params.random_seed:
            self._params.random_seed = random.randint(0, sys.maxint)
        random.seed(self._params.random_seed)
        if self._verbose:
            print "Random seed set to %d while building models" % self._params.random_seed

    def build_goniometer(self):

        self.goniometer = goniometer(matrix.col(self._params.goniometer.axis))

    def build_source(self):

        if self._params.beam.wavelength.random:
            wavelength = random.uniform(*self._params.beam.wavelength.range)
        else: wavelength = self._params.beam.wavelength.value

        assert self._params.beam.direction.method in ['inclination',
                                                      'close_to',
                                                      'exactly']

        if self._params.beam.direction.method == 'inclination':

            if self._params.beam.direction.inclination.random:
                inclination = random.gauss(0.,
                                  self._params.beam.direction.inclination.angle)
            else: inclination = self._params.beam.direction.inclination.angle

            beam = matrix.col((0, 0, 1)).rotate(matrix.col((0, 1, 0)),
                                                inclination, deg=True)

        elif self._params.beam.direction.method == 'close_to':

            beam = random_direction_close_to(
                matrix.col(self._params.beam.direction.close_to.direction),
                sd = self._params.beam.direction.close_to.sd)

        elif self._params.beam.direction.method == 'exactly':

            beam = matrix.col(self._params.beam.direction.exactly)

        self.source = source(beam, wavelength)

    def build_detector(self):

        assert self._params.detector.directions.method in ['close_to','exactly']

        if self._params.detector.directions.method == 'close_to':

            dir1 = random_direction_close_to(
                matrix.col(self._params.detector.directions.close_to.dir1),
                sd = self._params.detector.directions.close_to.sd)

            n = random_direction_close_to(
                matrix.col(self._params.detector.directions.close_to.norm),
                sd = self._params.detector.directions.close_to.sd)

        elif self._params.detector.directions.method == 'exactly':

            dir1 = matrix.col(self._params.detector.directions.exactly.dir1)

            n = matrix.col(self._params.detector.directions.exactly.norm)

        dir2 = n.cross(dir1).normalize()

        assert self._params.detector.centre.method in ['close_to','exactly']

        if self._params.detector.centre.method == 'close_to':

            centre = random_direction_close_to(
                matrix.col(self._params.detector.centre.close_to.value),
                sd = self._params.detector.centre.close_to.sd)

        elif self._params.detector.centre.method == 'exactly':

            centre = matrix.col(self._params.detector.centre.exactly.value)

        origin = centre - (0.5 * self._params.detector.npx_fast *
                           self._params.detector.pix_size * dir1 +
                           0.5 * self._params.detector.npx_slow *
                           self._params.detector.pix_size * dir2)
        self.detector = pilatus(origin, dir1, dir2,
                                    self._params.detector.pix_size,
                                    self._params.detector.pix_size,
                                    self._params.detector.npx_fast,
                                    self._params.detector.npx_slow)

    def _build_cell_vec(self, vec):

        if vec.length.random:
            length = random.uniform(*vec.length.range)
        else: length = vec.length.value

        assert vec.direction.method in ['close_to','exactly']

        if vec.direction.method == 'close_to':

            x = random_direction_close_to(
                    matrix.col(vec.direction.close_to.direction),
                    sd = vec.direction.close_to.sd)

        elif vec.direction.method == 'exactly':

            x = matrix.col(vec.direction.exactly.direction)

        return length * x



    def build_crystal(self):

        vecs = map(self._build_cell_vec,
                   [self._params.crystal.a,
                    self._params.crystal.b,
                    self._params.crystal.c])

        self.crystal = crystal(*vecs)
