#!/usr/bin/env python
#
# summation_3d_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface, Integration3DMixin

class Summation3dIntegrationExt(IntensityIface, Integration3DMixin):
  ''' Extension class to provide 3d summation integration. '''

  name = 'sum3d'

  def __init__(self, params, experiments):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration import Summation3d
    self._algorithm = Summation3d()

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

