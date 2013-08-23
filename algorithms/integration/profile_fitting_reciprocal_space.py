#
# profile_fitting_reciprocal_space.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.interfaces.integration import IntegrationInterface


class ProfileFittingReciprocalSpace(IntegrationInterface):
    '''A class to perform reciprocal space profile fitting'''

    def __init__(self, **kwargs):
        '''Initialise algorithm.'''

        # Set the parameters
        self.grid_size = kwargs['grid_size']
        self.threshold = kwargs['threshold']
        self.frame_interval = kwargs['frame_interval']
        self.bbox_nsigma = kwargs['n_sigma']

    def __call__(self, sweep, crystal, reflections, reference=None):
        '''Process the reflections.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflections to process

        Returns:
            The list of integrated reflections

        '''
        from dials.util.command_line import Command
        from dials.algorithms.reflection_basis import transform as rbt
        from dials.algorithms.integration import \
            ProfileFittingReciprocalSpaceAlgorithm

        # Initialise the reciprocal space transform
        Command.start('Initialising reciprocal space transform')
        transform = rbt.Forward(sweep, crystal, self.bbox_nsigma, self.grid_size)
        Command.end('Initialised reciprocal space transform')

        # Transform the reflections to reciprocal space
        Command.start('Transforming reflections to reciprocal space')
        transform(reflections)
        Command.end('Transformed {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Configure the integration algorithm with the locator class
        integrate = ProfileFittingReciprocalSpaceAlgorithm(reference)

        # Perform the integration
        Command.start('Integrating reflections in reciprocal space')
        integrate(reflections)
        Command.end('Integrated {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Return the integrated reflections
        return reflections
