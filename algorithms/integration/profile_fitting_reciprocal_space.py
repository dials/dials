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


class ProfileFittingReciprocalSpace(object):
  ''' Class to do reciprocal space profile fitting. '''

  def __init__(self, **kwargs):
    ''' Initialise the algorithm. '''
    from dials.algorithms.peak_finding.spot_matcher import SpotMatcher

    # Set the parameters
    self.grid_size = kwargs['grid_size']
    self.threshold = kwargs['threshold']
    self.frame_interval = kwargs['frame_interval']
    self.bbox_nsigma = kwargs['n_sigma']
    self.sigma_b = kwargs['sigma_b']
    self.sigma_m = kwargs['sigma_m']

  def __call__(self, experiment, reflections):
    ''' Do the integration.

    Transform the profiles to reciprocal space. Select the reflections
    to use in the reference profiles. Learn the reference profiles and
    integrate the intensities via profile fitting.

    '''
    assert("flags" in reflections)
    self._transform_profiles(experiment, reflections)
    self.learner = self._learn_references(experiment, reflections)
    return self._integrate_intensities(self.learner, reflections)

  def _transform_profiles(self, experiment, reflections):
    ''' Transform the reflection profiles to reciprocal space. '''
    from dials.util.command_line import Command
    from dials.algorithms.reflection_basis import transform as rbt
    from dials.array_family import flex

    # Initialise the reciprocal space transform
    Command.start('Initialising reciprocal space transform')
    spec = rbt.TransformSpec(experiment, self.sigma_b, self.sigma_m,
                             self.bbox_nsigma, self.grid_size)
    Command.end('Initialised reciprocal space transform')

    # Transform the reflections to reciprocal space
    Command.start('Transforming reflections to reciprocal space')
    s1 = reflections['s1']
    phi = reflections['xyzcal.mm'].parts()[2]
    shoebox = reflections['shoebox']
    rs_shoebox = flex.transformed_shoebox(spec, s1, phi, shoebox)
    reflections['rs_shoebox'] = rs_shoebox
    Command.end('Transformed {0} reflections'.format(len(reflections)))

  def _learn_references(self, experiment, reflections):
    ''' Learn the reference profiles. '''
    from dials.algorithms.integration.profile import GridSampler
    from dials.algorithms.integration.profile import GridReferenceLearner
    #from dials.algorithms.integration.profile import XdsCircleSampler
    #from dials.algorithms.integration.profile import XdsCircleReferenceLearner
    from dials.array_family import flex

    # Match the predictions with the strong spots
    #sind, pind = self.match(strong, reflections)
    pind = reflections.get_flags(reflections.flags.reference_spot)

    # Create the reference profile sampler
    assert(len(experiment.detector) == 1)
    image_size = experiment.detector[0].get_image_size()
    num_frames = experiment.scan.get_num_images()
    volume_size = image_size + (num_frames,)
    sampler = GridSampler(volume_size, (3, 3, 1))
    #sampler = XdsCircleSampler(volume_size, 1)

    # Configure the reference learner
    grid_size = (self.grid_size * 2 + 1,) * 3
    learner = GridReferenceLearner(sampler, grid_size, self.threshold)
    #learner = XdsCircleReferenceLearner(sampler, grid_size, self.threshold)
    profiles = reflections['rs_shoebox'].select(pind)
    coords = reflections['xyzcal.px'].select(pind)
    learner.learn(profiles, coords)

    # Return the learner
    return learner

  def _integrate_intensities(self, learner, reflections):
    ''' Integrate the intensities. '''
    from dials.util.command_line import Command
    from dials.algorithms.integration import ProfileFittingReciprocalSpaceAlgorithm

    # Configure the integration algorithm with the locator class
    integrate = ProfileFittingReciprocalSpaceAlgorithm(learner.locate())

    # Perform the integration
    Command.start('Integrating reflections in reciprocal space')
    profiles = reflections['rs_shoebox']
    coords = reflections['xyzcal.px']
    intensity = integrate(profiles, coords)
    mask = intensity.parts()[1] > 0
    reflections['intensity.raw.value'] = intensity.parts()[0]
    reflections['intensity.raw.variance'] = intensity.parts()[1]
    reflections['profile.correlation'] = intensity.parts()[2]
    reflections.set_flags(mask, reflections.flags.integrated)
    Command.end('Integrated {0} reflections'.format(len(reflections)))

    return reflections
