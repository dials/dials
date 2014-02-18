#!/usr/bin/env python
#
# kabsch_spotfinder_threshold.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces import SpotFinderThresholdIface


class KabschSpotFinderThresholdExt(SpotFinderThresholdIface):
  ''' Extensions to do xds-like threshold. '''

  name = 'xds'

  def __init__(self, params):
    ''' Initialise the algorithm. '''
    from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy
    self._algorithm = XDSThresholdStrategy(
      kernel_size=params.spotfinder.threshold.kernel_size,
      gain=params.lookup.gain_map,
      mask=params.lookup.mask,
      n_sigma_b=params.spotfinder.threshold.sigma_background,
      n_sigma_s=params.spotfinder.threshold.sigma_strong,
      min_count=params.spotfinder.threshold.min_local)

  def compute_threshold(self, image, mask):
    ''' Compute the threshold. '''
    return self._algorithm(image, mask)
