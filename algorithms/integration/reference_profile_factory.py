#!/usr/bin/env python
#
# dials.algorithms.integration.reference_profile_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class ReferenceProfileCreator(object):
    ''' A class to learn reference profiles from strong spots. '''

    def __init__(self, learn=None):
        ''' Initialise the algorithm.

        Params:
            learn The profile learning strategy

        '''
        from dials.algorithms.integration import ReflectionPredictor
        from dials.algorithms.peak_finding import SpotMatcher

        # Create the reflection predictor
        self.predict = ReflectionPredictor()

        # Create the spot matcher
        self.match = SpotMatcher()

        # Set the learning strategy
        self.learn = learn

    def __call__(self, sweep, crystal, strong):
        ''' Learn the reference profiles

        Params:
            sweep The sweep to process
            crystal The crystal to process
            strong The strong spots

        Returns:
            The reference profile locator

        '''
        # Predict the reflections
        predicted = self.predict(sweep, crystal)

        # Match the predictions with the strong spots
        matches = self.match(strong, predicted)

        # Learn the reference profiles
        return self.learn(sweep, crystal, matches)


class ProfileLearner(object):
    ''' A class to wrap the profile learning stuff. '''

    def __init__(self, bbox_nsigma, grid_size, threshold, frame_interval):
        ''' Configure the class.

        Params:
            bbox_nsigma The size of the bounding box
            grid_size The size of the transform grid size
            threshold The threshold level
            frame_interval The number of frames between references

        '''
        self.bbox_nsigma = bbox_nsigma
        self.frame_interval = frame_interval
        self.threshold = threshold
        self.grid_size = grid_size

    def __call__(self, sweep, crystal, reflections):
        ''' Learn the reference profiles

        Params:
            sweep The sweep class
            crystal The crystal class
            reflections The list of reflections to learn from

        Returns:
            The reference profile locator

        '''
        from dials.algorithms.reflection_basis import transform as rbt
        from dials.algorithms.integration.profile import XdsCircleSampler
        from dials.algorithms.integration.profile import XdsCircleReferenceLearner
        from dials.util.command_line import Command

        # Create the reference profile sampler
        image_size = sweep.get_detector().get_image_size()
        num_frames = sweep.get_scan().get_num_images()
        volume_size = image_size + (num_frames,)
        num_z = int(num_frames / self.frame_interval) + 1
        sampler = XdsCircleSampler(volume_size, num_z)

        # Configure the reference learner
        grid_size = (self.grid_size * 2 + 1,) * 3
        learner = XdsCircleReferenceLearner(sampler, grid_size, self.threshold)

        # Initialise the reciprocal space transform
        Command.start('Initialising reciprocal space transform')
        transform = rbt.Forward(sweep, crystal, self.bbox_nsigma, self.grid_size)
        Command.end('Initialised reciprocal space transform')

        # Transform the reflections to reciprocal space
        Command.start('Transforming reflections to reciprocal space')
        transform(reflections)
        Command.end('Transformed {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))

        # Learn the reference profiles
        Command.start('Learning reference profiles from {0} reflections'.format(
            len([r for r in reflections if r.is_valid()])))
        learner.learn(reflections)
        Command.end('Learnt {0} reference profiles from {1} reflections'.format(
            len(learner.locate()),
            len([r for r in reflections if r.is_valid()])))

        # Return the reference profile locator
        return learner.locate()


class ReferenceProfileFactory(object):
    ''' Factory class to create reference profile creators '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the integrator

        Params:
            params The input parameters

        Returns:
            The reference profile creator instance

        '''
        # Configure the strategies
        learn = ReferenceProfileFactory.configure_learner(params)

        # Return the reference profile creator with the given strategies
        return ReferenceProfileCreator(learn=learn)

    @staticmethod
    def configure_learner(params):
        ''' Configure the reference profile learner

        Params:
            params The input parameters

        Returns:
            The profile learner

        '''
        if params.integration.algorithm == "fit_rs":

            # Create the profile learner for reciprocal space fitting
            p = params.integration
            learner = ProfileLearner(
                bbox_nsigma=p.shoebox.n_sigma,
                grid_size=p.reciprocal_space.grid_size,
                threshold=p.profile.reference_signal_threshold,
                frame_interval=p.profile.reference_frame_interval)

        else:

            # If no profile fitting algorithm is selected, output a warning
            import sys
            import logging
            log = logging.getLogger(__name__)
            log.warn("No reference profiles were created for integration "
                     "algorithm '{0}'. In order to create reference profiles "
                     "you need to select a profile fitting integration "
                     "algorithm.".format(params.integration.algorithm))
            sys.exit(0)

        # Return the learner
        return learner
