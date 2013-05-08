#
# integrate3d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.interfaces.integration import IntegrationInterface

class Integrate3d(IntegrationInterface):
    '''A class to perform 3D integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        pass

    def __call__(self, reflections):
        '''Process the reflections.'''
        pass
