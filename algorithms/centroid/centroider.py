from __future__ import division
#
# centroider.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

class Centroider(object):
    ''' Class to centroid the reflections.'''

    def __init__(self):
        ''' Init the class '''
        pass

    def __call__(self, sweep, crystal, reflections):
        ''' Centroid the reflections.

        Params:
            sweep The sweep class
            crystal The crystal class
            reflections The reflection list

        Returns:
            The list of reflections

        '''
        from dials.algorithms.centroid import ComputeCentroid
        from dials.util.command_line import Command

        # Compute the reflections
        Command.start('Calculating reflection centroids')
        centroid = ComputeCentroid()
        centroid(reflections)
        Command.end('Calculated {0} reflection centroids'.format(
            len([r for r in reflections if r.is_valid()])))

        # Return the reflections
        return reflections
