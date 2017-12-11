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


class GaussianRSReferenceCalculatorFactory(object):
  '''
  A class to create the reference calculator

  '''

  @classmethod
  def create(Class,
             experiments,
             grid_size,
             scan_step,
             grid_method = 'circular_grid'):
    '''
    Create the intensity calculator

    '''
    from dials.algorithms.integration.parallel_integrator import GaussianRSReferenceCalculator
    from dials.algorithms.profile_model.modeller import SingleSampler
    from dials.algorithms.profile_model.modeller import CircleSampler
    from dials.algorithms.profile_model.modeller import GridSampler
    from dials.algorithms.profile_model.modeller import EwaldSphereSampler
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec

    from math import ceil

    # Assume the detector and scan are the same in each case
    detector = experiments[0].detector
    scan = experiments[0].scan

    # Get the number of scan points
    scan_range = scan.get_oscillation_range()
    scan_range = abs(scan_range[1] - scan_range[0])
    num_scan_points = int(ceil(scan_range / scan_step))

    # If multi panel then set to single
    if grid_method in ['regular_grid', 'circular_grid'] and len(detector) > 1:
      grid_method = 'single'

    # Create the sampler
    if grid_method == 'single':
      sampler = SingleSampler(
        scan.get_array_range(),
        num_scan_points)
    elif grid_method == 'regular_grid':
      sampler = GridSampler(
        detector[0].get_image_size(),
        scan.get_array_range(),
        (3, 3, num_scan_points))
    elif grid_method == 'circular_grid':
      sampler = CircleSampler(
        detector[0].get_image_size(),
        scan.get_array_range(),
        num_scan_points)
    elif grid_method == 'spherical_grid':
      sampler = EwaldSphereGridSampler(
        experiments[0].beam,
        experiments[0].detector,
        experiments[0].goniometer,
        experiments[0].scan,
        num_scan_points)
    else:
      raise RuntimeError("Unknown grid type")

    # Create the spec list
    spec_list = []
    for experiment in experiments:

      spec = TransformSpec(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        experiment.profile.sigma_b(deg=False),
        experiment.profile.sigma_m(deg=False),
        experiment.profile.n_sigma() * 1.5,
        grid_size)

      spec_list.append(spec)

    # Return the intensity algorithm
    return GaussianRSReferenceCalculator(sampler, spec_list)
