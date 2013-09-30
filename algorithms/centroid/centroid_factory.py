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
