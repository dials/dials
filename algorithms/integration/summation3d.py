#
# summation3d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.interfaces.integration import IntegrationInterface

class Summation3d(IntegrationInterface):
    '''A class to perform 3D summation integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''
        from dials.algorithms.integration import Summation3dAlgorithm
        self._integrator = Summation3dAlgorithm()

    def __call__(self, sweep, crystal, reflections):
        '''Process the reflections.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflections to process

        Returns:
            The list of integrated reflections

        '''
        from dials.util.command_line import Command

        # Integrate and return the reflections
        Command.start('Integrating reflections')
        self._integrator(reflections)
        Command.end('Integrated {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))
        return reflections
