#!/usr/bin/env python
#
# curved_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class CurvedBackgroundExt(BackgroundIface):
  ''' An extension class implementing curved background subtraction. '''

  name = 'curved'

  def __init__(self, params, experiment):
    ''' Initialise the algorithm. '''
    from dials.algorithms.background import CurvedSubtractor

    self._subtractor = CurvedSubtractor()

  def compute_background(self, reflections):
    ''' Compute the backgrond. '''
    from dials.util.command_line import Command

    # Do the background subtraction
    Command.start('Calculating reflection background')
    self._subtractor(reflections)
    Command.end('Calculated {0} background values'.format(len(reflections)))