from __future__ import absolute_import, division

import logging
logger = logging.getLogger("dials.extensions.global_spotfinder_threshold_ext")


class GlobalSpotFinderThresholdExt(object):
  ''' Extensions to do global threshold. '''

  name = 'single'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''
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
    self.params = params

  def compute_threshold(self, image, mask):
    '''
    Compute the threshold.

    :param image: The image to process
    :param mask: The pixel mask on the image
    :returns: A boolean mask showing foreground/background pixels

    '''
    return (image > self.params.spotfinder.threshold.single.global_threshold) & mask
