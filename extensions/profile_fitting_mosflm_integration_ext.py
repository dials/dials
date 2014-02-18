#!/usr/bin/env python
#
# profile_fitting_mosflm_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntegrationIface

class ProfileFittingMosflmIntegrationExt(IntegrationIface):
  ''' Extension class to provide mosflm profile fitting. '''

  name = 'mosflm'

  def __init__(self, params, experiment):
    ''' Initialise the algorithhm. '''
    from dials.algorithms.integration.mosflm_like import MosflmProfileFitting
    self._algorithm = MosflmProfileFitting(
        nblocks = params.integration.mosflm.nblocks)

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(reflections)
