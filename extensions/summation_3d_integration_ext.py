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

from dials.interfaces import IntegrationIface

class Summation3dIntegrationExt(IntegrationIface):
  ''' Extension class to provide 3d summation integration. '''

  name = 'sum3d'

  def __init__(self, params, reflections):
    ''' Initialise the algorithm. '''
    self._algorithm = Summation3d()

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

