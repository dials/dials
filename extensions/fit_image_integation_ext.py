#!/usr/bin/env python
#
# fit_image_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface

class FitImageIntegrationExt(IntensityIface):
  ''' Extension providing image space profile fitting. '''

  name = 'fit_image'

  phil = '''

    integrator = *auto 3d flat3d 2d single2d
      .type = choice
      .help = "The integrator to use."
      .expert_level=3

  '''

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration.fit_image import IntegrationAlgorithm
    self._algorithm = IntegrationAlgorithm(experiments, profile_model)

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

  @classmethod
  def type(cls, params, experiments):
    ''' Return the type of the integrator. '''
    from libtbx import Auto
    integrator_type = params.integration.intensity.sum.integrator
    if integrator_type == Auto or integrator_type == 'auto':
      integrator_type = '3d'
    return integrator_type
