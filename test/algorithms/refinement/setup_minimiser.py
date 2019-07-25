#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Setup experimental geometry for refinement test cases"""

# Python and cctbx imports
from __future__ import absolute_import, division, print_function
from libtbx.phil import parse, command_line

# Import the refinement engine
from dials.algorithms.refinement.engine import (
    SimpleLBFGS,
    LBFGScurvs,
    GaussNewtonIterations,
)


class Extract(object):
    """Parse and extract minimiser setup from PHIL"""

    def __init__(
        self,
        master_phil,
        target,
        prediction_parameterisation,
        local_overrides="",
        cmdline_args=None,
        verbose=True,
    ):

        self._target = target
        self._prediction_parameterisation = prediction_parameterisation
        self._verbose = verbose

        arg_interpreter = command_line.argument_interpreter(master_phil=master_phil)

        user_phil = parse(local_overrides)
        cmdline_phils = []
        if cmdline_args:
            for arg in cmdline_args:
                cmdline_phils.append(arg_interpreter.process(arg))

        working_phil = master_phil.fetch(sources=[user_phil] + cmdline_phils)

        self._params = working_phil.extract().minimiser.parameters

        self.refiner = self.build_minimiser()

    def build_minimiser(self):

        assert self._params.engine in ["SimpleLBFGS", "LBFGScurvs", "GaussNewton"]

        if self._params.engine == "SimpleLBFGS":
            refiner = SimpleLBFGS(
                target=self._target,
                prediction_parameterisation=self._prediction_parameterisation,
                log=self._params.logfile,
            )
            return refiner

        if self._params.engine == "LBFGScurvs":
            refiner = LBFGScurvs(
                target=self._target,
                prediction_parameterisation=self._prediction_parameterisation,
                log=self._params.logfile,
            )
            return refiner

        if self._params.engine == "GaussNewton":

            refiner = GaussNewtonIterations(
                target=self._target,
                prediction_parameterisation=self._prediction_parameterisation,
                log=self._params.logfile,
            )
            return refiner
