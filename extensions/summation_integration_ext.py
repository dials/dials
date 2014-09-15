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

from dials.interfaces import IntensityIface

class SummationIntegrationExt(IntensityIface):
  ''' Extension class to provide 3d summation integration. '''

  name = 'sum'

  phil = '''

    integrator = *3d flat3d 2d single2d
      .type = choice
      .help = "The integrator to use."
      .expert_level=3

  '''

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration import Summation
    self._algorithm = Summation()

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

  @classmethod
  def type(cls, params, experiments):
    ''' Return the type of the integrator. '''
    return params.integration.intensity.sum.integrator
