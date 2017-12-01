#!/usr/bin/env python
#
# algorithm.py
#
#  Copyright (C) 2017 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class GaussianRSMaskCalculatorFactory(object):
  '''
  Factory class for mask calculator

  '''

  @classmethod
  def create(Class, experiments):
    '''
    Create the mask calculator

    '''
    from dials.algorithms.integration.parallel_integrator import GaussianRSMaskCalculator
    from dials.algorithms.integration.parallel_integrator import GaussianRSMultiCrystalMaskCalculator
    result = GaussianRSMultiCrystalMaskCalculator()
    for e in experiments:
      alg = GaussianRSMaskCalculator(
        e.beam,
        e.detector,
        e.goniometer,
        e.scan,
        e.profile.delta_b(deg=False),
        e.profile.delta_m(deg=False))
      result.append(alg)
    return result


class GaussianRSIntensityCalculatorFactory(object):
  '''
  A class to create the intensity calculator

  '''

  @classmethod
  def create(Class,
             data,
             detector_space = False,
             deconvolution = False):
    '''
    Create the intensity calculator

    '''
    from dials.algorithms.integration.parallel_integrator import GaussianRSIntensityCalculator

    # Return the intensity algorithm
    return GaussianRSIntensityCalculator(data, detector_space, deconvolution)
