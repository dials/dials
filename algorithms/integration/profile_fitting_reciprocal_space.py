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

  reference_counter = 0

  def __init__(self, **kwargs):
    ''' Initialise the algorithm. '''
    from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
    from math import pi

    # Set the parameters
    self.grid_size = kwargs['grid_size']
    self.threshold = kwargs['threshold']
    self.frame_interval = kwargs['frame_interval']
    self.bbox_nsigma = kwargs['n_sigma']
    self.sigma_b = kwargs['sigma_b'] * pi / 180
    self.sigma_m = kwargs['sigma_m'] * pi / 180

  def __call__(self, experiments, reflections):
    ''' Do the integration.

    Transform the profiles to reciprocal space. Select the reflections
    to use in the reference profiles. Learn the reference profiles and
    integrate the intensities via profile fitting.

    '''
    from dials.model.serialize import dump
    assert(len(experiments) == 1)
    experiment = experiments[0]
    assert("flags" in reflections)
    assert(len(experiment.detector) == 1)
    self._integrate_by_summation(experiment, reflections)
    self._transform_profiles(experiment, reflections)
    self.learner = self._learn_references(experiment, reflections)
    counter = ProfileFittingReciprocalSpace.reference_counter
    dump.reference(self.learner.locate(), "reference_%d.pickle" % counter)
    ProfileFittingReciprocalSpace.reference_counter += 1
    return self._integrate_intensities(self.learner, reflections)

  def _integrate_by_summation(self, experiment, reflections):
    ''' Integrate the reflections by summation. '''
    from dials.algorithms.integration import Summation3d
    algorithm = Summation3d()
    algorithm(reflections)

  def _transform_profiles(self, experiment, reflections):
    ''' Transform the reflection profiles to reciprocal space. '''
    from dials.util.command_line import Command
    from dials.algorithms.reflection_basis import transform as rbt
    from dials.array_family import flex

    # Initialise the reciprocal space transform
    Command.start('Initialising reciprocal space transform')
    spec = rbt.TransformSpec(experiment, self.sigma_b, self.sigma_m,
                             5, self.grid_size)
    Command.end('Initialised reciprocal space transform')

    # Transform the reflections to reciprocal space
    Command.start('Transforming reflections to reciprocal space')
    s1 = reflections['s1']
    phi = reflections['xyzcal.mm'].parts()[2]
    shoebox = reflections['shoebox']
    rs_shoebox = flex.transformed_shoebox(spec, s1, phi, shoebox)
    reflections['rs_shoebox'] = rs_shoebox
    Command.end('Transformed {0} reflections'.format(len(reflections)))

    # Compute correlations between profile and idea
    Command.start('Compute correlation between profile and ideal')
    corr = rs_shoebox.correlation_with_ideal(5)
    reflections['correlation.ideal.profile'] = corr
    Command.end('Computed correlations between profile and ideal')

  def _learn_references(self, experiment, reflections):
    ''' Learn the reference profiles. '''
    from dials.algorithms.integration.profile import GridSampler
    from dials.algorithms.integration.profile import GridReferenceLearner
    from dials.algorithms.integration.profile import XdsCircleSampler
    from dials.algorithms.integration.profile import XdsCircleReferenceLearner
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
    counts = learner.counts()

    # Create the average profile
    from dials.util import pprint
    locator = learner.locate()

    print "Number of reflections contributing to each profile:"
    for i, num in enumerate(counts):
      xyz = locator.coord(i)
      print "%d: (%d, %d, %d); %d" % (
        i, int(xyz[0]), int(xyz[1]), int(xyz[2]), num)

    profiles = [locator.profile(i) for i in range(len(locator))]
    average_profile = sum(profiles) / len(profiles)
    print "Average Profile:\n"
    print pprint.profile3d(average_profile)

    # Return the learner
    return learner

  def _integrate_intensities(self, learner, reflections):
    ''' Integrate the intensities. '''
    from dials.array_family import flex
    from dials.util.command_line import Command
    from dials.algorithms.integration import ProfileFittingReciprocalSpaceAlgorithm

    # Configure the integration algorithm with the locator class
    integrate = ProfileFittingReciprocalSpaceAlgorithm(learner.locate())

    # Perform the integration
    Command.start('Integrating reflections in reciprocal space')
    mask = ~reflections.get_flags(reflections.flags.dont_integrate)
    profiles = reflections['rs_shoebox'].select(mask)
    coords = reflections['xyzcal.px'].select(mask)
    intensity = integrate(profiles, coords)
    I, I_var, P_cor = intensity.parts()
    mask2 = I_var > 0
    I_var.set_selected(mask2 != True, 0.0)
    reflections['intensity.prf.value'] = flex.double(len(reflections))
    reflections['intensity.prf.variance'] = flex.double(len(reflections))
    reflections['profile.correlation'] = flex.double(len(reflections))
    reflections['intensity.prf.value'].set_selected(mask, I)
    reflections['intensity.prf.variance'].set_selected(mask, I_var)
    reflections['profile.correlation'].set_selected(mask, P_cor)
    indices = flex.size_t(range(len(mask))).select(mask).select(~mask2)
    mask.set_selected(indices, False)
    reflections.set_flags(mask, reflections.flags.integrated_prf)
    Command.end('Integrated {0} reflections'.format(mask.count(True)))
    return reflections
