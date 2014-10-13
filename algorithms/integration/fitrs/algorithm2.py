#
# algorithm.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class IntegrationAlgorithm(object):
  ''' Class to do reciprocal space profile fitting. '''

  def __init__(self,
               experiments,
               profile_model,
               grid_size=5,
               threshold=2.0,
               debug=False):
    '''Initialise algorithm.'''
    assert(len(experiments) == len(profile_model))
    self._experiments = experiments
    self._profile_model = profile_model
    self._grid_size = grid_size
    self._threshold = threshold
    self._debug = debug

  def __call__(self, reflections):
    '''Process the reflections.

    Params:
        reflections The reflections to process

    Returns:
        The list of integrated reflections

    '''
    from dials.algorithms.integration.fitrs import ReciprocalSpaceProfileFitting
    from dials.algorithms.integration.fitrs import Spec
    from dials.algorithms.integration.interface import job_id
    from dials.array_family import flex
    from dials.util.command_line import Command

    # Get the flags
    flags = flex.reflection_table.flags

    # Create the algorithm
    algorithm = ReciprocalSpaceProfileFitting(self._grid_size)

    # Add the specs
    for experiment, model in zip(self._experiments, self._profile_model):
      algorithm.add(Spec(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        model.delta_b(deg=False),
        model.delta_m(deg=False)))

    # Perform the integration
    num = reflections.get_flags(flags.dont_integrate).count(False)
    Command.start('Integrating %d reflections with profile fitting' % num)
    profiles = algorithm.execute(reflections)

    mask1 = reflections.get_flags(reflections.flags.dont_integrate)
    mask2 = reflections.get_flags(reflections.flags.reference_spot)

    print mask1.count(True), mask2.count(True)

    # Print the number integrated
    num = reflections.get_flags(flags.integrated_prf).count(True)
    Command.end('Integrated %d reflections with profile fitting' % num)

    # Return the reflections
    return reflections
