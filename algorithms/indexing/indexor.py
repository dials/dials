#!/usr/bin/env python
#
# dials.algorithms.indexing.indexer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter, Nick Sauter, Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class IndexerSolution(object):
    '''FIXME need to define this so we can have a list of solutions from 
    indexing.'''

    pass

class Indexer(object):
    ''' The indexer base class. '''

    def __init__(self, strategies, parameters):
        ''' Initialise the indexer base class.

        Params:
            strategies: TBD; strategies for e.g. basis discovery
            parameters: TBD; the set of phil parameters

        '''
        # FIXME make these work properly by naming the strategies etc.
        self.strategies = strategies
        self.parameters = parameters
        # will work like this: 
        # self.do_this = do_this
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
        # - validate input (setup)
        # - discover beam centre (setup)
        # - map spots to RS (index)
        # - determine candidate basis vectors (index)
        # - determine basis sets ([P1_matrix], spots_and_lattice_id) (analyse)
        # - score possible lattice for each solution (analyse)
        # - refine lattice for each solution (refine)
        # - reject outliers for each solution (refine) 

        if while self._refined:
            while not self._analysed:
                while not self._indexed:
                    while not self._setuped:
                        self._setup()
                    self._index()
                self._analyse()
            self._refine()
        
        return

    def _setup(self):
        self._validate(self)
        self._discover_beam_centre_strategy()
        return

    def _index(self):
        self._index_strategy(self)
        return

    def _analyse(self):
        for lattice in self._lattices:
            self._analyse_strategy(lattice)
        return

    def _refine(self):
        for lattice in self._lattices:
            self._refine_strategy(lattice, spots)
        if not self._refined:
            # perhaps need to wind back to the mapping to reciprocal space and
            # try reindexing
            return
        for lattice in self._lattices:
            self._outlier_strategy(lattice, spots)
        # need to decide whether to reindexing the spot lists or what
        return
        
    def set_target_cell_lattice(self, cell, lattice):
        self._indexer_cell = cell
        self._indexer_lattice = lattice
        return

    def set_max_primitive_cell(self, max_primitive_cell):
        self._indexer_max_primitive_cell = max_primitive_cell
        return

    # etc.

class IndexerFactory(object):
    ''' Factory class to create indexers '''

    @staticmethod
    def get_from_somewhere_else(params):
        '''Get a different indexer implementation, which for example may 
        overload __call__ internally, or something'''

        # FIXME in here check with the registry for one of these based on the 
        # input PHIL parameters
        
        return False
    
    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the indexer

        Params:
            params The input parameters

        Returns:
            The indexer instance

        '''

        one_from_somewhere_else = IndexerFactory.get_from_somewhere_else(params)
        if one_from_somewhere_else:
            return one_from_somewhere_else

        # else configure a standard one with strategies
        
        # this is deliberately not implemented
        strategies = IndexerFactory.get_strategies_from_somewhere(params)
        
        # Return the indexer with the given strategies
        return Indexer(strategies, params)

    @staticmethod
    def get_strategies_from_somewhere(params):
        '''Get the strategies from somewhere, for example a registry.'''

        indexing = params.indexing

        return { } # or whatever; TBD

