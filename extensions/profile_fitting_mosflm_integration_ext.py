#!/usr/bin/env python
#
# profile_fitting_mosflm_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst & Luiso
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface, Integration3DMixin

class ProfileFittingMosflmIntegrationExt(IntensityIface, Integration3DMixin):
  ''' Extension class to provide mosflm profile fitting. '''

  name = 'mosflm'

  phil = '''
    nblocks = 5
      .help = "number of block per coordinate"
      .type = int
  '''

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithhm. '''
    from dials.algorithms.integration.mosflm_like import MosflmProfileFitting
    self._algorithm = MosflmProfileFitting(experiments,
        nblocks = params.integration.intensity.mosflm.nblocks)


  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)
