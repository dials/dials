#!/usr/bin/env cctbx.python

# Copyright (C) (2012) David Waterman, STFC Rutherford Appleton Laboratory, UK.
# This code is developed as part of the DIALS project and is provided for
# testing purposes only

"""Setup experimental geometry for refinement test cases"""

# Python and cctbx imports
from __future__ import division
from libtbx.phil import parse, command_line

# Import the refinement engine
from dials.scratch.dgw.refinement.engine import SimpleLBFGS, \
    LBFGScurvs, GaussNewtonIterations

class Extract(object):
    '''Parse and extract minimiser setup from PHIL'''

    def __init__(self, master_phil, target, prediction_parameterisation,
        local_overrides = "", cmdline_args = None, verbose=True):

        self._target = target
        self._prediction_parameterisation = prediction_parameterisation
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

        self._params = working_phil.extract().minimiser.parameters

        self.refiner = self.build_minimiser()

    def build_minimiser(self):

        assert self._params.engine in ["SimpleLBFGS", "LBFGScurvs",
            "GaussNewtonIterations"]

        if self._params.engine == "SimpleLBFGS":
            refiner = SimpleLBFGS(
                self._target,
                self._prediction_parameterisation,
                self._params.logfile,
                self._params.verbosity)
            return refiner

        if self._params.engine == "LBFGScurvs":
            refiner = LBFGScurvs(
                self._target,
                self._prediction_parameterisation,
                self._params.logfile,
                self._params.verbosity)
            return refiner

        if self._params.engine == "GaussNewtonIterations":

            refiner = GaussNewtonIterations(
                self._target,
                self._prediction_parameterisation,
                self._params.logfile,
                self._params.verbosity)
            return refiner
