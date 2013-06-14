#
# integrate2d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luiso & James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.interfaces.integration import IntegrationInterface

class Integrate2d(IntegrationInterface):
    '''A class to perform 2D integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        pass
    def __call__(self, sweep, crystal, reflections):
        '''Process the reflections.'''
        self.integrate(reflections)

    def integrate(self, reflections):
        from dials.algorithms.integration.summation2d \
         import tmp_translating_n_layering

        tmp_translating_n_layering(reflections)


