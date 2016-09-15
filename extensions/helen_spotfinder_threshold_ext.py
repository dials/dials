#!/usr/bin/env python
#
# helen_spotfinder_threshold.py
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
from logging import info


class HelenSpotFinderThresholdExt(SpotFinderThresholdIface):
  ''' Extensions to do spot finding threshold. '''

  name = 'helen'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''

      exp_spot_dimension = 3
        .type = int
        .help = "The expected spot dimensions"

      global_threshold = 100
        .type = float
        .help = "The global threshold value."

      min_blob_score = 0.7
        .type = float
        .help = "The minimum score for a blob"

      num_passes = 0
        .type = int
        .help = "Number of passes after updating ideal spot"

      debug = False
        .type = bool
        .help = "Write out correlation"

    ''')
    return phil

  def __init__(self, params):
    '''
    Initialise the algorithm.

    :param params: The input parameters

    '''
    self.params = params

  def compute_threshold(self, image, mask):
    '''
    Compute the threshold.

    :param image: The image to process
    :param mask: The pixel mask on the image
    :returns: A boolean mask showing foreground/background pixels

    '''
    from dials.algorithms.spot_finding.helen import BlobThresholdAlgorithm

    params = self.params.spotfinder.threshold.helen

    self._algorithm = BlobThresholdAlgorithm(
      pixels_per_row     = image.all()[1],
      row_count          = image.all()[0],
      exp_spot_dimension = params.exp_spot_dimension,
      global_threshold   = params.global_threshold,
      min_blob_score     = params.min_blob_score)
      num_passes         = params.num_passes)

    result = self._algorithm.threshold(image, mask)

    if self.params.spotfinder.threshold.helen.debug:
      from dials.array_family import flex
      corr = self._algorithm.correlation(image, mask)
      import cPickle as pickle
      pickle.dump(corr, open("correlation.pickle", "wb"))
    return result
