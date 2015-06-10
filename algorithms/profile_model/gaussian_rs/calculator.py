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
#
# FIXME Mosaicity seems to be overestimated
# FIXME Don't know how XDS REFLECTING_RANGE is calculated
# FIXME Don't know what XDS REFLECTION_RANGE is used for
# FIXME Don't know what XDS BEAM_DIVERGENCE is used for
# FIXME Should maybe be scan varying
# FIXME Don't know how XDS calculated the n_sigma

from __future__ import division


class ComputeEsdBeamDivergence(object):
  '''Calculate the E.s.d of the beam divergence.'''

  def __init__(self, detector, reflections):
    ''' Calculate the E.s.d of the beam divergence.

    Params:
        detector The detector class
        reflections The reflections

    '''
    from scitbx.array_family import flex
    from math import sqrt

    # Calculate the beam direction variances
    variance = self._beam_direction_variance_list(detector, reflections)

    # Calculate and return the e.s.d of the beam divergence
    self._sigma = sqrt(flex.sum(variance) / len(variance))

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

    # Get the reflection columns
    shoebox = reflections['shoebox']
    bbox = reflections['bbox']
    xyz = reflections['xyzobs.px.value']

    # Loop through all the reflections
    variance = []
    for r in range(len(reflections)):

      # Get the coordinates and values of valid shoebox pixels
      mask = shoebox[r].mask != 0
      coords = shoebox[r].coords(mask)
      values = shoebox[r].values(mask)
      s1 = shoebox[r].beam_vectors(detector, mask)

      # Calculate the beam vector at the centroid
      panel = shoebox[r].panel
      s1_centroid = detector[panel].get_pixel_lab_coord(xyz[r][0:2])
      angles = s1.angle(s1_centroid, deg=False)
      variance.append(flex.sum(values * (angles**2)) / (flex.sum(values) - 1))

    # Return a list of variances
    return flex.double(variance)


class FractionOfObservedIntensity(object):
  '''Calculate the fraction of observed intensity for different sigma_m.'''

  def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
    '''Initialise the algorithm. Calculate the list of tau and zetas.

    Params:
        reflections The list of reflections
        experiment The experiment object

    '''
    from dials.array_family import flex
    from math import sqrt

    # Get the oscillation width
    dphi2 = scan.get_oscillation(deg=False)[1] / 2.0

    # Calculate a list of angles and zeta's
    tau, zeta = self._calculate_tau_and_zeta(crystal, beam, detector,
                                             goniometer, scan, reflections)

    # Calculate zeta * (tau +- dphi / 2) / sqrt(2)
    self.e1 = (tau + dphi2) * flex.abs(zeta) / sqrt(2.0)
    self.e2 = (tau - dphi2) * flex.abs(zeta) / sqrt(2.0)

  def _calculate_tau_and_zeta(self, crystal, beam, detector, goniometer, scan, reflections):
    '''Calculate the list of tau and zeta needed for the calculation.

    Params:
        reflections The list of reflections
        experiment The experiment object.

    Returns:
        (list of tau, list of zeta)

    '''
    from scitbx.array_family import flex

    # Calculate the list of frames and z coords
    bbox = reflections['bbox']
    phi = reflections['xyzcal.mm'].parts()[2]

    # Calculate the zeta list
    zeta = reflections['zeta']

    # Calculate the list of tau values
    tau = []
    zeta2 = []
    scan = scan
    for b, p, z in zip(bbox, phi, zeta):
      for f in range(b[4], b[5]):
        phi0 = scan.get_angle_from_array_index(int(f), deg=False)
        phi1 = scan.get_angle_from_array_index(int(f)+1, deg=False)
        tau.append((phi1 + phi0) / 2.0 - p)
        zeta2.append(z)

    # Return the list of tau and zeta
    return flex.double(tau), flex.double(zeta2)

  def __call__(self, sigma_m):
    '''Calculate the fraction of observed intensity for each observation.

    Params:
        sigma_m The mosaicity

    Returns:
        A list of log intensity fractions

    '''
    from math import sqrt
    from scitbx.array_family import flex
    import scitbx.math

    # Tiny value
    TINY = 1e-10
    assert(sigma_m > TINY)

    # Calculate the two components to the fraction
    a = scitbx.math.erf(self.e1 / sigma_m)
    b = scitbx.math.erf(self.e2 / sigma_m)

    # Calculate the fraction of observed reflection intensity
    R = (a - b) / 2.0

    # Set any points <= 0 to 1e-10 (otherwise will get a floating
    # point error in log calculation below).
    assert(R.all_ge(0))
    mask = R < TINY
    assert(mask.count(True) < len(mask))
    R.set_selected(mask, TINY)

    # Return the logarithm of r
    return flex.log(R)


class ComputeEsdReflectingRange(object):
  '''calculate the e.s.d of the reflecting range (mosaicity).'''


  class Estimator(object):
    '''Estimate E.s.d reflecting range by maximum likelihood estimation.'''

    def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
      '''Initialise the optmization.'''
      from scitbx import simplex
      from scitbx.array_family import flex
      from math import pi, exp
      import random

      # Initialise the function used in likelihood estimation.
      self._R = FractionOfObservedIntensity(crystal, beam, detector, goniometer,
                                            scan, reflections)

      # Set the starting values to try
      start = random.random() * pi / 180
      stop = random.random() * pi / 180
      starting_simplex = [flex.double([start]), flex.double([stop])]

      # Initialise the optimizer
      optimizer = simplex.simplex_opt(
        1,
        matrix=starting_simplex,
        evaluator=self,
        tolerance=1e-7)

      # Get the solution
      self.sigma = exp(optimizer.get_solution()[0])

    def target(self, log_sigma):
      ''' The target for minimization. '''
      from scitbx.array_family import flex
      from math import exp
      return -flex.sum(self._R(exp(log_sigma[0])))

  def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
    '''initialise the algorithm with the scan.

    params:
        scan the scan object

    '''

    # Calculate sigma_m
    estimator = ComputeEsdReflectingRange.Estimator(
      crystal, beam, detector, goniometer, scan, reflections)

    # Save the solution
    self._sigma = estimator.sigma

  def sigma(self):
    ''' Return the E.S.D reflecting rang. '''
    return self._sigma


class ProfileModelCalculator(object):
  ''' Class to help calculate the profile model. '''

  def __init__(self, reflections, crystal, beam, detector, goniometer, scan, min_zeta=0.05):
    ''' Calculate the profile model. '''
    from logging import info
    from dxtbx.model.experiment.experiment_list import Experiment
    from dials.array_family import flex

    # Check input has what we want
    assert(reflections is not None)
    assert("miller_index" in reflections)
    assert("s1" in reflections)
    assert("shoebox" in reflections)
    assert("xyzobs.px.value" in reflections)
    assert("xyzcal.mm" in reflections)

    # Calculate the E.S.D of the beam divergence
    info('Calculating E.S.D Beam Divergence.')
    beam_divergence = ComputeEsdBeamDivergence(detector, reflections)

    # Set the sigma b
    self._sigma_b = beam_divergence.sigma()

    # FIXME Calculate properly
    if goniometer is None or scan is None:
      self._sigma_m = 0.00001
    else:

      # Select by zeta
      zeta = reflections.compute_zeta(
        Experiment(
          crystal=crystal,
          beam=beam,
          detector=detector,
          goniometer=goniometer,
          scan=scan))
      mask = flex.abs(zeta) >= min_zeta
      reflections = reflections.select(mask)

      # Calculate the E.S.D of the reflecting range
      info('Calculating E.S.D Reflecting Range.')
      reflecting_range = ComputeEsdReflectingRange(crystal, beam, detector,
                                                   goniometer, scan, reflections)

      # Set the sigmas
      self._sigma_m = reflecting_range.sigma()

  def sigma_b(self):
    ''' Return the E.S.D beam divergence. '''
    return self._sigma_b

  def sigma_m(self):
    ''' Return the E.S.D reflecting range. '''
    return self._sigma_m

class ScanVaryingProfileModelCalculator(object):
  ''' Class to help calculate the profile model. '''

  def __init__(self, reflections, crystal, beam, detector, goniometer, scan, min_zeta=0.05):
    ''' Calculate the profile model. '''
    from logging import info
    from copy import deepcopy
    from collections import defaultdict
    from dials.array_family import flex
    from dxtbx.model.experiment.experiment_list import Experiment

    # Check input has what we want
    assert(reflections is not None)
    assert("miller_index" in reflections)
    assert("s1" in reflections)
    assert("shoebox" in reflections)
    assert("xyzobs.px.value" in reflections)
    assert("xyzcal.mm" in reflections)

    # Select by zeta
    zeta = reflections.compute_zeta(
      Experiment(
        crystal=crystal,
        beam=beam,
        detector=detector,
        goniometer=goniometer,
        scan=scan))
    mask = flex.abs(zeta) >= min_zeta
    reflections = reflections.select(mask)

    # Split the reflections into partials
    reflections = deepcopy(reflections)
    reflections.split_partials_with_shoebox()

    # Get a list of reflections for each frame
    bbox = reflections['bbox']
    index_list = defaultdict(list)
    for i, (x0, x1, y0, y1, z0, z1) in enumerate(bbox):
      assert(z1 == z0 + 1)
      index_list[z0].append(i)
    reflection_list = dict((key, reflections.select(flex.size_t(value)))
                           for key, value in index_list.iteritems())

    # The range of frames
    z0, z1 = scan.get_array_range()
    min_z = min(index_list.iterkeys())
    max_z = max(index_list.iterkeys())
    assert(z0 == min_z)
    assert(z1 == max_z+1)

    # Compute for all frames
    self._num = []
    sigma_b = []
    sigma_m = []
    for i in range(z0, z1):

      # Get reflections at the index
      reflections = reflection_list[i]

      self._num.append(len(reflections))

      info('Computing profile model for frame %d' % i)

      # Calculate the E.S.D of the beam divergence
      beam_divergence = ComputeEsdBeamDivergence(detector, reflections)

      # Set the sigma b
      sigma_b.append(beam_divergence.sigma())

      # Calculate the E.S.D of the reflecting range
      reflecting_range = ComputeEsdReflectingRange(crystal, beam, detector,
                                                   goniometer, scan, reflections)

      # Set the sigmas
      sigma_m.append(reflecting_range.sigma())

    def convolve(data, kernel):
      assert(len(kernel) & 1)
      result = []
      mid = len(kernel) // 2
      for i in range(len(data)):
        r = 0
        for j in range(len(kernel)):
          k = i - mid + j
          if k < 0:
            r += kernel[j] * data[0]
          elif k >= len(data):
            r += kernel[j] * data[-1]
          else:
            r += kernel[j] * data[k]
        result.append(r)
      return result

    def gaussian_kernel(n):
      from math import exp
      assert(n & 1)
      mid = n // 2
      sigma = mid / 3.0
      kernel = []
      for i in range(n):
        kernel.append(exp(-(i-mid)**2 / (2*sigma**2)))
      kernel = [k / sum(kernel) for k in kernel]
      return kernel

    # Smooth the parameters
    kernel = gaussian_kernel(51)
    sigma_b_new = convolve(sigma_b, kernel)
    sigma_m_new = convolve(sigma_m, kernel)
    self._sigma_b = flex.double(sigma_b_new)
    self._sigma_m = flex.double(sigma_m_new)
    assert(len(self._sigma_b) == len(self._sigma_m))

  def num(self):
    ''' The number of reflections used. '''
    return self._num

  def sigma_b(self):
    ''' Return the E.S.D beam divergence. '''
    return self._sigma_b

  def sigma_m(self):
    ''' Return the E.S.D reflecting range. '''
    return self._sigma_m
