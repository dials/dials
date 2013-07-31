#
# masker.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class Masker(object):
    '''A class to perform all the shoebox masking.'''

    def __init__(self, sweep, detector_mask, delta_d, delta_m):
        ''' Initialise the masking algorithms

        Params:
            sweep The sweep object
            detector_mask The detector mask
            delta_d The extent of the reflection in reciprocal space
            delta_m The extent of the reflection in reciprocal space

        '''
        from dials.algorithms.shoebox import MaskBadPixels
        from dials.algorithms.shoebox import MaskOverlapping
        from dials.algorithms.shoebox import MaskForeground

        # Construct the bad pixel mask
        self.mask_bad_pixels = MaskBadPixels(detector_mask)

        # Construct the overlapping reflection mask
        self.mask_overlapping = MaskOverlapping()

        # Construct the foreground pixel mask
        self.mask_foreground = MaskForeground(
            sweep.get_beam(), sweep.get_detector(),
            sweep.get_goniometer(), sweep.get_scan(),
            delta_d, delta_m)

    def __call__(self, reflections, adjacency_list=None):
        ''' Mask the given reflections.

        Params:
            reflections The reflection list
            adjacency_list The adjacency_list (optional)

        Returns:
            The masked reflection list

        '''
        from dials.util.command_line import Command

        # Mask the bad pixels
        Command.start('Masking bad pixels')
        self.mask_bad_pixels(reflections)
        Command.end('Masked bad pixels for {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Mask the overlaps if an adjacency list is given
        if adjacency_list:
            Command.start('Masking overlapping reflections')
            self.mask_overlapping(reflections, adjacency_list)
            Command.end('Masked {0} overlapping reflections'.format(
                len(adjacency_list)))

        # Mask the foreground region
        Command.start('Masking foreground pixels')
        self.mask_foreground(reflections)
        Command.end('Masked foreground pixels for {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Return the reflections
        return reflections
