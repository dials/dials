#
# summation_reciprocal_space.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.interfaces.integration import IntegrationInterface

class SummationReciprocalSpace(IntegrationInterface):
    '''A class to perform 3D summation integration'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''

        # Get the parameters
        self._bbox_nsigma = kwargs['n_sigma']
        self._grid_size = kwargs['grid_size']

    def __call__(self, sweep, crystal, reflections):
        '''Process the reflections.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflections to process

        Returns:
            The list of integrated reflections

        '''
        from dials.algorithms.reflection_basis import transform as rbt
        from dials.algorithms.integration import \
            SummationReciprocalSpaceAlgorithm
        from dials.util.command_line import Command

        # Initialise the reciprocal space transform
        Command.start('Initialising reciprocal space transform')
        transform = rbt.Forward(
            sweep.get_beam(),
            sweep.get_detector(),
            sweep.get_goniometer(),
            sweep.get_scan(),
            crystal.get_mosaicity(),
            self._bbox_nsigma,
            self._grid_size)
        Command.end('Initialised reciprocal space transform')

        # Initialise the integration algorithm
        Command.start('Initialising integration algorithm')
        integrate = SummationReciprocalSpaceAlgorithm(
            sweep.get_beam(),
            sweep.get_detector(),
            sweep.get_goniometer())
        Command.end('Initialied integration algorithm')

        # Transform the reflections to reciprocal space
        Command.start('Transforming reflections to reciprocal space')
        transform(reflections)
        Command.end('Transformed {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Integrate the reflections
        Command.start('Integrate the reflections in reciprocal space')
        integrate(reflections)
        Command.end('Integrated {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Return the reflections
        return reflections
