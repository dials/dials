#
# summation3d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Summation3d(object):
  '''A class to perform 3D summation integration'''

  def __init__(self, **kwargs):
    '''Initialise algorithm.'''
    pass

  def __call__(self, experiment, reflections, reference=None):
    '''Process the reflections.

    Params:
        experiment The experiment data
        reflections The reflections to process

    Returns:
        The list of integrated reflections

    '''
    from dials.algorithms.integration import integrate_by_summation
    from dials.util.command_line import Command

    # Integrate and return the reflections
    Command.start('Integrating reflections')
    intensity = reflections['shoebox'].summed_intensity_foreground()
    reflections['intensity.raw.value'] = intensity.observed_value()
    reflections['intensity.raw.variance'] = intensity.observed_variance()
    Command.end('Integrated {0} reflections'.format(len(reflections)))
    return reflections
