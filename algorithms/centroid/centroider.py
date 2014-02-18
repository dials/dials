#
# centroider.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Centroider(object):
  ''' Class to centroid the reflections.'''

  def __init__(self):
    ''' Init the class '''
    pass

  def __call__(self, experiment, reflections):
    ''' Centroid the reflections.

    Params:
        experiment The experiment data
        reflections The reflection list

    Returns:
        The list of reflections

    '''
    from dials.util.command_line import Command

    # Compute the reflections
    Command.start('Calculating reflection centroids')
    reflections['shoebox'].centroid_valid()
    Command.end('Calculated {0} reflection centroids'.format(len(reflections)))

    # Return the reflections
    return reflections
