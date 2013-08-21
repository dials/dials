#!/usr/bin/env python
#
# dials.algorithms.integration.reference_profile_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class ReferenceProfileCreator(object):
    ''' A class to learn reference profiles from strong spots. '''

    def __init__(self):
        ''' Initialise the algorithm.'''
        from dials.algorithms.integration import ReflectionPredictor
        from dials.algorithms.peak_finding import SpotMatcher

        # Create the reflection predictor
        self.predict = ReflectionPredictor()

        # Create the spot matcher
        self.match = SpotMatcher()

        self.learn = lambda x: []

    def __call__(self, sweep, crystal, strong):
        ''' Learn the reference profiles

        Params:
            sweep The sweep to process
            crystal The crystal to process
            strong The strong spots

        '''

        # Predict the reflections
        predicted = self.predict(sweep, crystal)

        # Match the predictions with the strong spots
        matches = self.match(strong, predicted)

        # Learn the reference profiles
        return self.learn(matches)


class ReferenceProfileFactory(object):
    ''' Factory class to create reference profile creators '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the integrator

        Params:
            params The input parameters

        Returns:
            The reference profile creator instance

        '''
        # Return the reference profile creator with the given strategies
        return ReferenceProfileCreator()
