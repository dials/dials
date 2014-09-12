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

  phil = '''
    kernel_size = 3 3
      .help = "The size of the local area around the spot in which"
              "to calculate the mean and variance. The kernel is"
              "given as a box of size (2 * nx + 1, 2 * ny + 1) centred"
              "at the pixel."
      .type = ints(size=2)
      .expert_level = 1

    sigma_background = 6
      .help = "The number of standard deviations of the coefficient of"
              "variation (variance / mean) in the local area below"
              "which the pixel will be classified as background."
      .type = float
      .expert_level = 1

    sigma_strong = 3
      .help = "The number of standard deviations above the mean in the"
              "local area above which the pixel will be classified as"
              "strong."
      .type = float
      .expert_level = 1

    min_local=2
      .help = "The number of pixels in the local area of each pixel needed"
              "to do the thresholding. Setting to 0 or less means that all"
              "the pixels under the kernel are needed. The minimum allowable"
              "number is 2"
      .type = int
      .expert_level = 1
  '''

  default = True

  def __init__(self, params):
    ''' Initialise the algorithm. '''
    from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy
    self._algorithm = XDSThresholdStrategy(
      kernel_size=params.spotfinder.threshold.xds.kernel_size,
      gain=params.spotfinder.lookup.gain_map,
      mask=params.spotfinder.lookup.mask,
      n_sigma_b=params.spotfinder.threshold.xds.sigma_background,
      n_sigma_s=params.spotfinder.threshold.xds.sigma_strong,
      min_count=params.spotfinder.threshold.xds.min_local)

  def compute_threshold(self, image, mask):
    ''' Compute the threshold. '''
    return self._algorithm(image, mask)
