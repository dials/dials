#!/usr/bin/env python
#
# profile_model.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ComputeEsdBeamDivergence(object):
  '''Calculate the E.s.d of the beam divergence.'''

  def __init__(self, detector, reflections):
    ''' Calculate the E.s.d of the beam divergence.

    Params:
        detector The detector class
        reflections The reflections

    '''
    # Calculate the beam direction variances
    variance = self._beam_direction_variance_list(detector, reflections)

    # Calculate and return the e.s.d of the beam divergence
    self._sigma = self._compute_esd(variance)

  def sigma(self):
    ''' Return the E.S.D of the beam divergence. '''
    return self._sigma

  def _beam_direction_variance_list(self, detector, reflections):
    '''Calculate the variance in beam direction for each spot.

    Params:
        reflections The list of reflections

    Returns:
        The list of variances

    '''
    from scitbx.array_family import flex
    from scitbx import matrix
    import numpy

    # Get the reflection columns
    shoebox = reflections['shoebox']
    bbox = reflections['bbox']
    xyz = reflections['xyzobs.px.value']

    # Loop through all the reflections
    variance = []
    for r in range(len(reflections)):

      # Find the active pixels
      zi, yi, xi = numpy.where(shoebox[r].mask.as_numpy_array() != 0)
      index = zip(map(int, zi), map(int, yi), map(int, xi))

      # Extract the pixel values
      values = flex.double([shoebox[r].data[k, j, i] for k, j, i in index])

      # Get the pixel coordinates centres
      xs, xf, ys, yf, zs, zf = bbox[r]
      xp = xi + xs + 0.5
      yp = yi + ys + 0.5
      zp = zi + zs + 0.5

      # Calculate the beam directions to each pixel
      s1 = [detector[0].get_pixel_lab_coord((x, y)) for x, y in zip(xp, yp)]

      # Calculate the beam vector at the centroid
      xc, yc, zc = xyz[r]
      s1_centroid = detector[0].get_pixel_lab_coord((xc, yc))

      # Calculate the variance in beam vector directions
      var = self._beam_direction_variance(s1_centroid, s1, values)
      variance.append(var)

    # Return a list of variances
    return variance

  def _beam_direction_variance(self, s1_centroid, s1, values):
    '''Calculate the angles between the s1 centroid and s1 vectors for
    each pixel and then calculate the variance of the angles.

    Params:
        s1_centroid The centroid beam direction
        s1 The list of beam directions
        values The list of pixel values

    Returns:
        The variance in the beam direction

    '''
    from scitbx.array_family import flex
    from scitbx import matrix

    # Calculate angles between vectors
    s1_centroid = matrix.col(s1_centroid)
    angles = flex.double([s1_centroid.angle(matrix.col(s), deg=True) for s in s1])

    # Calculate the variance of the angles
    return flex.sum(values * (angles**2)) / flex.sum(values)

  def _compute_esd(self, variance):
    '''Calculate the beam divergence as the sum of centroid variance of the
    intensity weighted diffracted beam directionsl.

    Params:
        variance The variance of the beam directions

    Returns:
        The e.s.d of the beam divergence

    '''
    from math import sqrt

    # Return the beam divergence as the sum / num reflections
    return sqrt(sum(variance) / len(variance))


class ProfileModel(object):
  ''' Class to help calculate the profile model. '''

  def __init__(self, experiment, reflections):
    ''' Calculate the profile model. '''

    # Check input has what we want
    assert(experiment is not None)
    assert(reflections is not None)
    assert("miller_index" in reflections)
    assert("s1" in reflections)
    assert("shoebox" in reflections)
    assert("xyzobs.px.value" in reflections)

    # Calculate the E.S.D of the beam divergence
    beam_divergence = ComputeEsdBeamDivergence(experiment.detector, reflections)


    self._sigma_b = beam_divergence.sigma()
    self._sigma_m = 0

  def sigma_b(self):
    return self._sigma_b

  def sigma_m(self):
    return self._sigma_m
