#!/usr/bin/env python
#
# fitrs_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface

class FitrsIntegrationExt(IntensityIface):
  ''' Extension providing reciprocal space profile fitting. '''

  name = 'fitrs'

  default = True

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''

      grid_size = 5
        .type = int
        .help = "The size of the profile grid."

      threshold = 0.02
        .type = float
        .help = "The threshold to use in reference profile"

      single_reference = False
        .type = bool
        .help = "Use a single reference profile for profile fitting"

      debug = False
        .type = bool
        .help = "Save the reference profiles and other debug info."

    ''')
    return phil

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration.fitrs.algorithm import IntegrationAlgorithm
    self._algorithm = IntegrationAlgorithm(
      experiments,
      profile_model,
      grid_size=params.integration.intensity.fitrs.grid_size,
      threshold=params.integration.intensity.fitrs.threshold,
      single_reference=params.integration.intensity.fitrs.single_reference,
      debug=params.integration.intensity.fitrs.debug)

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)

  @classmethod
  def type(cls, params, experiments):
    ''' Return the type of the integrator. '''
    return '3d'

