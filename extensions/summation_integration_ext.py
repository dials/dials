#!/usr/bin/env python
#
# summation_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface, Integration3DMixin

class SummationIntegrationExt(IntensityIface, Integration3DMixin):
  ''' Extension class to provide 3d summation integration. '''

  name = 'sum'

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration import Summation
    self._algorithm = Summation()

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

