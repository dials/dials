#!/usr/bin/env python
#
# dials.process.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class PredictAndMatch(object):
    ''' A class to predict and match with observations. '''

    def __init__(self):
        ''' Initialise the algorithm.

        Params:
            learn The profile learning strategy

        '''
        from dials.algorithms.integration import ReflectionPredictor
        from dials.algorithms.peak_finding import SpotMatcher

        # Create the reflection predictor
        self.predict = ReflectionPredictor()

        # Create the spot matcher
        self.match = SpotMatcher()

    def __call__(self, sweep, crystal, observed):
        ''' Learn the reference profiles

        Params:
            sweep The sweep to process
            crystal The crystal to process
            observed The strong spots

        Returns:
            The matched reflections

        '''
        # Predict the reflections
        predicted = self.predict(sweep, crystal)

        # Match the predictions with the strong spots
        return self.match(observed, predicted)


class Process(object):
    ''' The processing pipeline '''

    def __init__(self, find_spots, refine_geometry, create_reference, integrate):
        ''' Initialise the dials pipeline with the required strategies

        Params:
            find_spots Find the strong spots
            refine_geometry Refine the geometry
            create_reference Create the reference profiles
            integrate Integrate the reflection intensities

        '''
        # Set all the strategies
        self.find_spots_internal = find_spots
        self.refine_geometry_internal = refine_geometry
        self.create_reference_internal = create_reference
        self.integrate_internal = integrate

        # Set all the properties
        self.observed = None
        self.predicted = None
        self.reference = None
        self.integrated = None

    def find_spots(self, sweep):
        ''' A wrapper around the find spots function '''
        print '\n' + '=' * 80
        print 'Find Spots'
        print '=' * 80 + '\n'
        return self.find_spots_internal(sweep)

    def refine_geometry(self, sweep, crystal):
        ''' A wrapper around the refine geometry function '''
        print '\n' + '=' * 80
        print 'Refine Geometry'
        print '=' * 80 + '\n'

        # Do the initial prediction and match with observations
        match_observations = PredictAndMatch()
        observed = match_observations(sweep, crystal, self.observed)

        # Do the finement
        print "\nPerforming Refinement\n"
        self.refine_geometry_internal(sweep, crystal, observed)

        # Predict some new spots
        return self.refine_geometry_internal.predict_reflections()

    def create_reference(self, sweep, crystal):
        ''' A wrapper around the refernece profiles function '''
        print '\n' + '=' * 80
        print 'Creating reference profiles'
        print '=' * 80 + '\n'
        return self.create_reference_internal(sweep, crystal, self.observed)

    def integrate(self, sweep, crystal):
        ''' A wrapper around the integrate function '''
        print '\n' + '=' * 80
        print 'Integrating Reflections'
        print '=' * 80 + '\n'
        return self.integrate_internal(sweep, crystal, reference=self.reference)

    def run(self, sweep, crystal):
        ''' Perform all the processing.

        Params:
            sweep The sweep to process
            crystal The crystal to process

        Returns:
            The integrated intensities

        '''
        # Process the sweep and get the integrated intensities
        # Do this in the following way:
        #  - First of all, find the strong spots.
        #  - Perform the indexing
        #  - Refine the experimental geometry
        #  - Create the reference profiles
        #  - Integrated the reflections
        self.observed   = self.find_spots(sweep)
        self.predicted  = self.refine_geometry(sweep, crystal)
        self.reference  = self.create_reference(sweep, crystal)
        self.integrated = self.integrate(sweep, crystal)


class ProcessFactory(object):
    ''' A factory to create the processing pipeline '''

    @staticmethod
    def from_parameters(params, verbosity):
        ''' Construct the pipeline from the input parameters

        Params:
            params The input parameters
            verbosity Verbosity for refinement

        Returns:
            The processing pipeline

        '''
        from dials.algorithms.peak_finding.spotfinder_factory \
            import SpotFinderFactory as SpotFinderFactory
        from dials.algorithms.refinement import RefinerFactory
        from dials.algorithms.integration import ReferenceProfileFactory
        from dials.algorithms.integration import IntegratorFactory

        # Get the spot finder from the input parameters
        print 'Configuring spot finder from input parameters'
        find_spots = SpotFinderFactory.from_parameters(params)

        # Get the refiner from the input parameters
        print 'Configuring refiner from input parameters'
        refine_geometry = RefinerFactory.from_parameters(params, verbosity)

        # Get the reference profile creator from the input parameters
        print 'Configuring reference profile creator from input parameters'
        create_reference = ReferenceProfileFactory.from_parameters(params)

        # Get the integrator from the input parameters
        print 'Configurating integrator from input parameters'
        integrate = IntegratorFactory.from_parameters(params)

        # Construct and return the pipeline
        return Process(find_spots=find_spots,
                       refine_geometry=refine_geometry,
                       create_reference=create_reference,
                       integrate=integrate)
