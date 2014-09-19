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
  ''' A class to perform profile fitting '''

  def __init__(self, experiments, **kwargs):
    '''Initialise algorithm.'''
    self._experiments = experiments

  def __call__(self, reflections):
    '''Process the reflections.

    Params:
        reflections The reflections to process

    Returns:
        The list of integrated reflections

    '''
    from dials.algorithms.integration.fit_image import ImageSpaceProfileFitting
    from dials.algorithms.integration.fit_image import Spec
    from dials.array_family import flex

    # Get the flags
    flags = flex.reflection_table.flags

    # Create the algorithm
    algorithm = ImageSpaceProfileFitting()

    # Add the specs
    for experiment in self._experiments:
      algorithm.add(Spec(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan))

    # Perform the integration
    algorithm.execute(reflections)

    # Print the number integrated
    num = reflections.get_flags(flags.integrated_prf).count(True)
    print 'Integated %d reflections with profile fitting' % num

    # Return the reflections
    return reflections
