#
# summation.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class IntegrationAlgorithm(object):
  '''A class to perform 3D summation integration'''

  def __init__(self, **kwargs):
    '''Initialise algorithm.'''
    pass

  def __call__(self, reflections):
    '''Process the reflections.

    Params:
        experiment The experiment data
        reflections The reflections to process

    Returns:
        The list of integrated reflections

    '''
    from logging import info
    from dials.array_family import flex
    from time import time

    # Integrate and return the reflections
    start_time = time()
    info('')
    info(' Beginning integration by summation')
    info('  using %d reflections' % len(reflections))
    intensity = reflections['shoebox'].summed_intensity()
    reflections['intensity.sum.value'] = intensity.observed_value()
    reflections['intensity.sum.variance'] = intensity.observed_variance()
    success = intensity.observed_success()
    reflections.set_flags(success, reflections.flags.integrated_sum)
    info('  successfully processed %d reflections' % success.count(True))
    info('  time taken: %g seconds' % (time() - start_time))
    return reflections
