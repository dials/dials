#!/usr/bin/env python
#
# flat_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class FlatBackgroundExt(BackgroundIface):
  ''' An extension class implementing flat background subtraction. '''

  name = 'flat'

  def __init__(self, params, experiment):
    ''' Initialise the algorithm. '''
    from dials.algorithms.background import FlatSubtractor

    self._subtractor = FlatSubtractor()

  def compute_background(self, reflections):
    ''' Compute the backgrond. '''
    from dials.util.command_line import Command

    # Do the background subtraction
    #next line temporarily commented for consistency (need to consult James)
    #Command.start('Calculating reflection background')
    self._subtractor(reflections)
    #next line temporarily commented for consistency (need to consult James)
    #Command.end('Calculated {0} background values'.format(len(reflections)))
