#!/usr/bin/env python
#
# gain.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

class ComputeBackgroundAndGain(object):
  '''Class to calculate the gain map of a detector.'''

  def __init__(self, mask, kernel_size=(3, 3)):
    '''Initialise the class with a mask of trusted regions'''
    from scitbx.array_family import flex
    self._mask = mask.as_1d().as_int()
    self._mask.reshape(flex.grid(mask.all()))
    self._gain = None
    self._background = None
    self._kernel_size = kernel_size
    self._count = 1

  def add(self, image):
    '''Add an image to calculation.'''

    from scitbx.array_family import flex
    import numpy

    # Get the mask as an integer array
    mask = self._mask.deep_copy()

    # Calculate the gain map
    gain_map, mask, mean = self._compute_gain_map(image, mask)

    # Get the mask with strong pixels
    mask = self._remove_strong_pixels(gain_map, mask)

    # Set all masked indices to 1
    gain_map = gain_map.as_numpy_array()
    index = numpy.where(mask.as_numpy_array() == 0)
    gain_map[index] = 1.0
    gain_map = flex.double(gain_map)

    # Add the gain to the sum
    if self._gain == None:
      self._gain = gain_map
    else:
      self._gain += gain_map

    # Add the image to the background sum
    if self._background == None:
      self._background = mean
    else:
      self._background += mean

    # Add to the count
    self._count += 1

  def _compute_gain_map(self, image, mask):
    '''Complete the gain map of a single image.'''
    from dials.algorithms.image.filter import index_of_dispersion_filter

    # Need double image
    image = image.as_double()

    # Filter the image and return the mask
    filteralg = index_of_dispersion_filter(image, mask, self._kernel_size, 0)
    filtered  = filteralg.index_of_dispersion()
    mask      = filteralg.mask()
    mean      = filteralg.mean()

    # Return the filtered image and mask
    return filtered, mask, mean

  def _remove_strong_pixels(self, gain_map, mask):
    '''Remove strong pixels from the gain map.'''
    from scitbx.array_family import flex
    from operator import mul
    from math import sqrt
    import numpy

    # Elements in kernel
    n = reduce(mul, self._kernel_size)

    # Get only those pixels that are poisson distributed
    bound = 1.0 + 3.0 * (2.0 / sqrt(n - 1))
    index = numpy.where(gain_map.as_numpy_array() > bound)
    mask = mask.as_numpy_array()
    mask[index] = 0
    mask = flex.int(mask)

    # Return the mask
    return mask

  def gain(self):
    '''Compute the full gain map.'''
    from dials.algorithms.image.filter import mean_filter

    # Divide all gain values by count and smooth the gain
    self._gain /= self._count
    self._gain *= self._mask.as_double()
    return mean_filter(self._gain, self._mask, self._kernel_size, 0)

  def background(self):
    '''Compute the full background map.'''
    from dials.algorithms.image.filter import mean_filter
    self._background /= self._count
    self._gain *= self._mask.as_double()
    return mean_filter(self._background, self._mask,
        self._kernel_size, 0)
