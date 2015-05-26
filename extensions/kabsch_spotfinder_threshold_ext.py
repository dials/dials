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
from libtbx.utils import Sorry
from dials.interfaces import SpotFinderThresholdIface


class KabschSpotFinderThresholdExt(SpotFinderThresholdIface):
  ''' Extensions to do xds-like threshold. '''

  name = 'xds'

  default = True

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''
      gain = None
        .type = float(value_min=0.0)
        .help = "Use a flat gain map for the entire detector. Cannot be used"
                "in conjunction with lookup.gain_map parameter."

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

      min_local = 2
        .help = "The minimum number of pixels under the image processing kernel"
                "that are need to do the thresholding operation. Setting the"
                "value between 2 and the total number of pixels under the"
                "kernel will force the algorithm to use that number as the"
                "minimum. If the value is less than or equal to zero, then"
                "the algorithm will use all pixels under the kernel. In"
                "effect this will add a border of pixels which are always"
                "classed as background around the edge of the image and around"
                "any masked out pixels."
        .type = int
        .expert_level = 1

      global_threshold = 0
        .type = float
        .help = "The global threshold value. Consider all pixels less than this"
                "value to be part of the background."
    ''')
    return phil

  def __init__(self, params):
    '''
    Initialise the algorithm.

    :param params: The input parameters

    '''
    from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy
    self._algorithm = XDSThresholdStrategy(
      kernel_size=params.spotfinder.threshold.xds.kernel_size,
      gain=params.spotfinder.threshold.xds.gain,
      mask=params.spotfinder.lookup.mask,
      n_sigma_b=params.spotfinder.threshold.xds.sigma_background,
      n_sigma_s=params.spotfinder.threshold.xds.sigma_strong,
      min_count=params.spotfinder.threshold.xds.min_local,
      global_threshold=params.spotfinder.threshold.xds.global_threshold)

  def compute_threshold(self, image, mask):
    '''
    Compute the threshold.

    :param image: The image to process
    :param mask: The pixel mask on the image
    :returns: A boolean mask showing foreground/background pixels

    '''
    return self._algorithm(image, mask)
