#
# __init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class IntegrationInterface(object):
    '''Interface for integration algorithms.'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        pass

    def __call__(self, reflections):
        '''Process the reflections.'''
        raise RuntimeError('Overload!')
