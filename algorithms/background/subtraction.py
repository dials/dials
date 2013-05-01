#
# dials.algorithms.background.subtraction.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces.background import BackgroundSubtractionInterface

class BackgroundSubtractor(BackgroundSubtractionInterface):
    '''A class to perform background subtraction.'''

    def __init__(self, **kwargs):
        '''Initialise the algorithm.'''

        from dials.algorithms.integration import SubtractBackground

        # Create the object to do the background subtraction
        self._subtract_background = SubtractBackground(
            kwargs.get('min_pixels', 10),
            kwargs.get('num_sigma', 3))

    def __call__(self, reflections):
        '''Subtract the background the reflections.'''
        self._subtract_background(reflections)
