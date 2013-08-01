#!/usr/bin/env python
#
# dials.algorithms.integration.integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class CentroidRunner(object):
    ''' The centroid runner class. '''

    def __init__(self, compute_spots, compute_centroid):
        ''' Initialise the integrator base class.

        Params:
            compute_spots The spot extractor strategy
            compute_centroid The centroid strategy

        '''
        self.compute_spots = compute_spots
        self.compute_centroid = compute_centroid

    def __call__(self, sweep, crystal, predicted=None):
        ''' Call to calculate centroids.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            predicted The predicted reflection list

        Returns:
            A reflection list

        '''
        # Extract the reflections from the sweep
        if predicted == None:
            predicted = self.compute_spots(sweep, crystal)

        # Calculate the reflection centroids
        return self.compute_centroid(sweep, crystal, predicted)


class CentroidFactory(object):
    ''' Factory class to create centroiders '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the centroid algorithm

        Params:
            params The input parameters

        Returns:
            The integrator instance

        '''
        # Configure the algorithms to predict reflections and calculate
        # the centroids needed for refinement etc.
        compute_spots = CentroidFactory.configure_predictor(params)
        compute_centroid = CentroidFactory.configure_centroid(params)

        # Return the centroider with the given strategies
        return CentroidRunner(compute_spots=compute_spots,
                              compute_centroid=compute_centroid)

    @staticmethod
    def configure_predictor(params):
        ''' Given a set of parameters, configure the reflection predictor

        Params:
            params The input parameters

        Returns:
            The extractor instance

        '''
        from dials.algorithms.integration import ReflectionPredictor

        # Return the reflection extractor instance
        return ReflectionPredictor()

    @staticmethod
    def configure_centroid(params):
        ''' Given a set of parameters, configure the centroid calculator

        Params:
            params The input parameters

        Returns:
            The centroid calculator instance

        '''
        from dials.algorithms.peak_finding.strong_spot_picker \
            import StrongSpotOrganiser

        # Load some lookup maps
        gain_map = CentroidFactory.load_image(params.lookup.gain_map)
        dark_map = CentroidFactory.load_image(params.lookup.dark_map)

        # Create the threshold strategy
        threshold = CentroidFactory.configure_threshold(params, gain_map)

        # Create the centroider
        return StrongSpotOrganiser(
            min_spot_size=params.spotfinder.filter.min_spot_size,
            max_separation=params.spotfinder.filter.max_separation,
            threshold_strategy=threshold)

    @staticmethod
    def configure_threshold(params, gain_map):
        '''Get the threshold strategy'''
        from dials.algorithms.peak_finding.threshold \
            import UnimodalThresholdStrategy, XDSThresholdStrategy

        # Chose the strategy
        if params.spotfinder.threshold.algorithm == 'xds':
            return XDSThresholdStrategy(
                kernel_size=params.spotfinder.threshold.kernel_size,
                gain=gain_map,
                n_sigma_b=params.spotfinder.threshold.sigma_background,
                n_sigma_s=params.spotfinder.threshold.sigma_strong)

        elif params.spotfinder.threshold.algorithm == 'unimodal':
            return UnimodalThresholdStrategy()

        else:
            raise RuntimeError('Unknown threshold strategy')

    @staticmethod
    def load_image(filename):
        ''' Given a filename, load an image

        Params:
            filename The input filename

        Returns:
            The image or None

        '''
        from dials.util import image

        # If no filename is set then return None
        if not filename:
            return None

        # Read the image and return the image data
        handle = image.reader()
        handle.read_file(filename)
        return handle.get_data()
