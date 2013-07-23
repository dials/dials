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
        self.n_div = kwargs['n_div']


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
        from dials.algorithms.reflection_basis import transform as rbt
        from dials.algorithms.integration.profile import \
            XdsCircleSampler
        from dials.algorithms.integration.profile import \
            XdsCircleReferenceLearner
        from dials.algorithms.integration import \
            ProfileFittingReciprocalSpaceAlgorithm

        # Create the reference profile sampler
        image_size = sweep.get_detector().get_image_size()
        num_frames = sweep.get_scan().get_num_images()
        volume_size = image_size + (num_frames,)
        num_z = int(num_frames / self.frame_interval) + 1
        sampler = XdsCircleSampler(volume_size, num_z)

        # Initialise the reciprocal space transform
        Command.start('Initialising reciprocal space transform')
        transform = rbt.Forward(
            sweep.get_beam(),
            sweep.get_detector(),
            sweep.get_scan(),
            crystal.get_mosaicity(),
            self.bbox_nsigma,
            self.grid_size)
        Command.end('Initialised reciprocal space transform')

        # Transform the reflections to reciprocal space
        Command.start('Transforming reflections to reciprocal space')
        transform(reflections)
        Command.end('Transformed {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Configure the reference learner
        grid_size = (self.grid_size * 2 + 1,) * 3
        learner = XdsCircleReferenceLearner(
            sampler, grid_size, self.threshold)

        # Learn the reference profiles around the detector
        Command.start('Learning reference profiles')
        learner.learn(reflections)
        locate = learner.locate()
        Command.end('Learnt {0} reference profiles'.format(locate.size()))

        # Configure the integration algorithm with the locator class
        integrate = ProfileFittingReciprocalSpaceAlgorithm(locate)

        # Perform the integration
        Command.start('Integrating reflections in reciprocal space')
        integrate(reflections)
        Command.end('Integrated {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Return the integrated reflections
        return reflections
