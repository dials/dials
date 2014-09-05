#!/usr/bin/env python
#
# summation_2d_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface, Integration3DMixin

class Summation2dIntegrationExt(IntensityIface, Integration3DMixin):
  ''' Extension to provide 2d summation integration. '''

  name = 'sum2d'

  def __init__(self, params, experiment):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration import Summation2d
    self._algorithm = Summation2d()

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

