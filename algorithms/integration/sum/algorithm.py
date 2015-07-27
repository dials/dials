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
    '''
    Initialise algorithm.

    '''
    pass

  def __call__(self, reflections):
    '''Process the reflections.

    :param reflections: The reflections to integrate
    :return: The list of integrated reflections

    '''
    from dials.array_family import flex

    # Integrate and return the reflections
    intensity = reflections['shoebox'].summed_intensity()
    reflections['intensity.sum.value'] = intensity.observed_value()
    reflections['intensity.sum.variance'] = intensity.observed_variance()
    success = intensity.observed_success()
    reflections.set_flags(success, reflections.flags.integrated_sum)
    return success
