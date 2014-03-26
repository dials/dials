#!/usr/bin/env python
#
# xds_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class XdsBackgroundExt(BackgroundIface):
  ''' An extension class implementing XDS background subtraction. '''

  name = 'xds'

  phil = '''
    min_pixels = 10
      .help = "The minimum number of pixels to use in calculating the"
              "background intensity."
      .type = int
  '''

  def __init__(self, params, experiment):
    ''' Initialise the algorithm. '''
    from dials.algorithms.background import XdsSubtractorAlgorithm

    if params:
      min_data = params.integration.background.xds.min_pixels
    else:
      min_data = 10

    self._subtractor = XdsSubtractorAlgorithm(min_data)

  def compute_background(self, reflections):
    ''' Compute the backgrond. '''
    from dials.util.command_line import Command

    # Do the background subtraction
    Command.start('Calculating reflection background')
    mask = self._subtractor(reflections['shoebox'])
    reflections.del_selected(mask != True)
    Command.end('Calculated {0} background values'.format(len(reflections)))
