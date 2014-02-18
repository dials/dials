#
# threshold.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ThresholdStrategy(object):
  '''Base class for spot finder threshold strategies.'''

  def __init__(self, **kwargs):
    '''Initialise with key word arguments.'''
    pass

  def __call__(self, image):
    '''Threshold the image.'''
    raise RuntimeError('Overload Me!')


class UnimodalThresholdStrategy(ThresholdStrategy):
  '''Unimodal histogram thresholding strategy.'''

  def __init__(self, **kwargs):
    '''Initialise the threshold.'''

    # Initialise the base class
    ThresholdStrategy.__init__(self, **kwargs)

    # Get the arguments
    trusted_range = kwargs.get('trusted_range', (0, 20000))

    # Make sure the range is valid
    self._hrange = (0, int(trusted_range[1]))

  def __call__(self, image):
    '''Calculate the threshold for this image.

    Params:
        image The image to process
        trusted_range The trusted range of pixel values

    Returns:
        The threshold value

    '''
    from dials.algorithms.image.threshold import maximum_deviation
    from dials.algorithms.image.threshold import probability_distribution

    # Get the probability distribution from the image
    p = probability_distribution(image, self._hrange)

    # Calculate the threshold and add to list
    threshold = maximum_deviation(p)

    # Return a threshold mask
    return image >= threshold


class XDSThresholdStrategy(ThresholdStrategy):
  '''A class implementing a 'gain' threshold.'''

  def __init__(self, **kwargs):
    '''Set the threshold algorithm up

    Params:
        kwargs The keyword arguments

    '''
    # Initialise the base class
    ThresholdStrategy.__init__(self, **kwargs)

    # Get the parameters
    self._kernel_size = kwargs.get('kernel_size', (3, 3))
    self._gain        = kwargs.get('gain')
    self._n_sigma_b   = kwargs.get('n_sigma_b', 6)
    self._n_sigma_s   = kwargs.get('n_sigma_s', 3)
    self._min_count   = kwargs.get('min_count', 2)

  def __call__(self, image, mask):
    '''Call the thresholding function

    Params:
        image The image to process

    Returns:
        The threshold image mask

    '''
    from dials.algorithms.image import threshold

    # Do the thresholding, if gain is given then use gain threshold,
    # otherwise do normal poisson exclusion (fano) threshold
    result = None
    if self._gain:
      result = threshold.kabsch_w_gain(
        image.as_double(),
        mask,
        self._gain,
        self._kernel_size,
        self._n_sigma_b,
        self._n_sigma_s,
        self._min_count)
    else:
      result = threshold.kabsch(
        image.as_double(),
        mask,
        self._kernel_size,
        self._n_sigma_b,
        self._n_sigma_s,
        self._min_count)

    # Return the result
    return result
