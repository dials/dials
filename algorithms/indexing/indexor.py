#!/usr/bin/env python
#
# dials.algorithms.indexing.integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter, Nick Sauter, Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Indexor(object):
    ''' The indexor base class. '''

    def __init__(self, strategies, parameters):
        ''' Initialise the indexor base class.

        Params:
            strategies: TBD; strategies for e.g. basis discovery
            parameters: TBD; the set of phil parameters

        '''
        # FIXME make these work properly by naming the strategies etc.
        self.strategies = strategies
        self.parameters = parameters
        # will work like this: 
        # self.correct_intensity = correct_intensity
        # however N.B. that this will probably mean you cannot pickle an
        # instance of this

    def __call__(self, spots, detectors, beam, goniometer = None, scan = None):
        ''' Call to index.

        Params:
            spots: The spot list inc. detector number
            detectors: dxtbx detector list
            beam: beam information
            goniometer: the goniometer; optional (for e.g. still images)
            scan: the scan information; optional (for e.g. still images)

        Returns:
            TBD

        '''

        # structured loops within loops to employ input strategies to - 
        # 
        # - validate input
        # - discover beam centre
        # - map spots to RS
        # - determine candidate basis vectors
        # - determine basis sets => ([P1_A_matrix], indexed_spots_and_lattice_id)
        # - score possible lattice for each solution
        # - refine lattice for each solution
        # - reject outliers for each solution

class IndexorFactory(object):
    ''' Factory class to create indexors '''

    @staticmethod
    def get_from_somewhere_else(params):
        '''Get a different indexor implementation, which for example may 
        overload __call__ internally, or something'''

        # FIXME in here check with the registry for one of these based on the 
        # input PHIL parameters
        
        return False
    
    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the indexor

        Params:
            params The input parameters

        Returns:
            The indexor instance

        '''

        one_from_somewhere_else = IndexorFactory.get_from_somewhere_else(params)
        if one_from_somewhere_else:
            return one_from_somewhere_else

        # else configure a standard one with strategies
        
        # this is deliberately not implemented
        strategies = IndexorFactory.get_strategies_from_somewhere(params)
        
        # Return the indexor with the given strategies
        return Indexor(strategies, params)

    @staticmethod
    def get_strategies_from_somewhere(params):
        '''Get the strategies from somewhere, for example a registry.'''

        indexing = params.indexing

        return { } # or whatever; TBD

