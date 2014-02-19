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
    return flex.sum(values * (angles**2)) / (flex.sum(values) - 1)

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


class FractionOfObservedIntensity(object):
  '''Calculate the fraction of observed intensity for different sigma_m.'''

  def __init__(self, reflections, experiment):
    '''Initialise the algorithm. Calculate the list of tau and zetas.

    Params:
        reflections The list of reflections
        experiment The experiment object

    '''
    self.dphi = experiment.scan.get_oscillation(deg=False)[1]
    self.tau, self.zeta = self._calculate_tau_and_zeta(reflections, experiment)

  def _calculate_tau_and_zeta(self, reflections, experiment):
    '''Calculate the list of tau and zeta needed for the calculation.

    Params:
        reflections The list of reflections
        experiment The experiment object.

    Returns:
        (list of tau, list of zeta)

    '''
    from scitbx.array_family import flex
    from dials.algorithms.reflection_basis import zeta_factor

    # Calculate the list of frames and z coords
    bbox = reflections['bbox']
    frames = [range(b[4], b[5]) for b in bbox]
    phi = reflections['xyzcal.mm'].parts()[2]
    s1s = reflections['s1']

    # Calculate the zeta list
    s0 = experiment.beam.get_s0()
    m2 = experiment.goniometer.get_rotation_axis()
    zeta = [zeta_factor(s0, s1, m2) for s1 in s1s]

    # Calculate the list of tau values
    tau = []
    zeta2 = []
    scan = experiment.scan
    for rf, p, z in zip(frames, phi, zeta):
      for f in rf:
        phi0 = scan.get_angle_from_array_index(int(f), deg=False)
        phi1 = scan.get_angle_from_array_index(int(f)+1, deg=False)
        tau.append((phi1 + phi0) / 2.0 - p)
        zeta2.append(z)

    # Return the list of tau and zeta
    assert(len(zeta2) == len(tau))
    return flex.double(tau), flex.double(zeta2)

  def __call__(self, sigma_m):
    '''Calculate the fraction of observed intensity for each observation.

    Params:
        sigma_m The mosaicity

    Returns:
        A list of fractions of length n

    '''
    from math import sqrt, erf, exp
    from scitbx.array_family import flex
    import numpy

    # Tiny value
    TINY = 1e-10

    # Oscillation range / 2
    dphi2 = self.dphi / 2

    # Calculate the denominator to the fraction
    den =  sqrt(2) * exp(sigma_m) / flex.abs(self.zeta)

    # Calculate the two components to the fraction
    a = flex.double([erf(x) for x in (self.tau + dphi2) / den])
    b = flex.double([erf(x) for x in (self.tau - dphi2) / den])

    # Calculate the fraction of observed reflection intensity
    R = (a - b) / 2.0

    # Set any points <= 0 to 1e-10 (otherwise will get a floating
    # point error in log calculation below).
    bad_index = numpy.where(R.as_numpy_array() < TINY)[0]
    for i in bad_index:
      R[int(i)] = TINY

    # Return the logarithm of r
    return flex.log(R)


class MaximumLikelihoodEstimator(object):
  '''Estimate E.s.d reflecting range by maximum likelihood estimation.'''

  def __init__(self, reflections, experiment):
    '''Initialise the optmization.'''
    from scitbx import simplex
    from scitbx.array_family import flex
    from math import pi
    import random

    # Initialise the function used in likelihood estimation.
    self._R = FractionOfObservedIntensity(reflections, experiment)

    # Set the starting values to try
    start = random.random() * pi / 180
    stop = random.random() * pi / 180
    starting_simplex = [flex.double([start]), flex.double([stop])]

    # Initialise the optimizer
    self._optimizer = simplex.simplex_opt(1, matrix=starting_simplex,
        evaluator=self, tolerance=1e-7)

  def target(self, sigma):
    '''Return the target function evaluated at this valid.

    Params:
        sigma The parameters

    Returns:
        The value of the target function at the given value.

    '''
    return -self.likelihood(sigma[0])

  def likelihood(self, sigma_m):
    '''Return the likelihood of the current sigma

    Params:
        sigma_m the estimated mosaicity

    Returns:
        The likelihood of this value.

    '''
    from scitbx.array_family import flex
    return flex.sum(self._R(sigma_m))

  def solve(self):
    '''Perform maximum likelihood estimation of sigma_m

    Returns:
        The value of sigma_m

    '''
    from math import exp
    return exp(self._optimizer.get_solution()[0])



class ComputeEsdReflectingRange(object):
  '''calculate the e.s.d of the reflecting range (mosaicity).'''

  def __init__(self, experiment, reflections):
    '''initialise the algorithm with the scan.

    params:
        scan the scan object

    '''
    maximum_likelihood = MaximumLikelihoodEstimator(
        reflections, experiment)
    self._sigma = maximum_likelihood.solve()

  def sigma(self):
    ''' Return the E.S.D reflecting rang. '''
    return self._sigma


class ProfileModel(object):
  ''' Class to help calculate the profile model. '''

  def __init__(self, experiment, reflections):
    ''' Calculate the profile model. '''
    from dials.util.command_line import Command
    from math import pi

    # Check input has what we want
    assert(experiment is not None)
    assert(reflections is not None)
    assert("miller_index" in reflections)
    assert("s1" in reflections)
    assert("shoebox" in reflections)
    assert("xyzobs.px.value" in reflections)

    # Calculate the E.S.D of the beam divergence
    Command.start('Calculating E.S.D Beam Divergence.')
    beam_divergence = ComputeEsdBeamDivergence(experiment.detector, reflections)
    Command.end('Calculated E.S.D Beam Divergence')

    # Calculate the E.S.D of the reflecting range
    Command.start('Calculating E.S.D Reflecting Range.')
    reflecting_range = ComputeEsdReflectingRange(experiment, reflections)
    Command.end('Calculated E.S.D Reflecting Range.')

    # Set the sigmas
    self._sigma_b = beam_divergence.sigma()
    self._sigma_m = reflecting_range.sigma() * 180.0 / pi

  def sigma_b(self):
    ''' Return the E.S.D beam divergence. '''
    return self._sigma_b

  def sigma_m(self):
    ''' Return the E.S.D reflecting range. '''
    return self._sigma_m
