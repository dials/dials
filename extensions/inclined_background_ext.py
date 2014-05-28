#!/usr/bin/env python
#
# inclined_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class InclinedBackgroundExt(BackgroundIface):
  ''' An extension class implementing inclined background subtraction. '''

  name = 'inclined'

  def __init__(self, params, experiment):
    ''' Initialise the algorithm. '''
    from dials.algorithms.background import InclinedSubtractor

    self._subtractor = InclinedSubtractor()

  def compute_background(self, reflections):
    ''' Compute the backgrond. '''
    from dials.util.command_line import Command

    # Do the background subtraction

    #next line temporarily commented for consistency (need to consult James)
    #Command.start('Calculating reflection background')
    self._subtractor(reflections)
    #next line temporarily commented for consistency (need to consult James)
    #Command.end('Calculated {0} background values'.format(len(reflections)))
