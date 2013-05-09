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
        return self._subtract_background(reflections)


class PixelDiscrimination(BackgroundSubtractionInterface):
    '''Class to perform pixel discrimination background subtraction.

    This class performs background subtraction by first selecting which
    pixels are peak and which are background. It then calculates the
    background and subtracts it from the reflection pixels. Both methods
    are set as strategies in the constructor with keyword arguments
      - discriminate - The pixel discrimination strategy
      - subtract - The background subtraction strategy

    '''

    def __init__(self, **kwargs):
        '''Initialise the algorithm.

        Params:
            kwargs - The keyword arguments

        '''
        # Get the pixel discrimination and subtraction stratgies
        self._discriminate = kwargs.get('discriminate', NormalDiscriminator())
        self._subtract     = kwargs.get('subtract', MeanSubtractor())

    def __call__(self, reflections):
        '''Perform background subtraction on all reflections

        Params:
            reflections The list of reflections

        Return:
            Array of booleans True/False for success/no success

        '''
        from scitbx.array_family import flex

        # For each reflection assign the pixels in the shoebox area
        # to be either peak or background. In input, any non-zero mask
        # pixels are said to belong to the reflection. On output
        # -       0  mask pixels don't belong to the reflection
        # - (1 << 0) mask pixels are background
        # - (1 << 1) mask pixels are peak
        # Indeterminate pixels can optionally be marked both peak
        # and background
        success1 = self._discriminate(reflections)

        # Get indices and reduced list of reflections
        index = flex.select(range(len(success1)), flags=success1)
        reflections2 = ReflectionList(flex.select(reflections, flags=success1))

        # Having calculated which pixels belong to the peak and which
        # pixels belong to the background. The background is calculated
        # and then subtracted from the peak pixels.
        success2 = self._subtract(reflections2)

        # Update the arrays of success/fails
        for i, idx in enumerate(index):
            success1[idx] = success2[i]

        # Return whether it was successful
        return success1
