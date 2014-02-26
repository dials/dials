#!/usr/bin/env python
#
# lui_spotfinder_threshold.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuente-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces import SpotFinderThresholdIface


class SmoothingSpotFinderThresholdExt(SpotFinderThresholdIface):
  ''' Extensions to do luiso's' threshold. '''

  name = 'Smoothing'

  def __init__(self, params):
    ''' Initialise the algorithm. '''
    pass

  def compute_threshold(self, image, mask):
    ''' Compute the threshold. '''
    raise RuntimeError('Not Implemented')
