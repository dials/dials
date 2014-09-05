#!/usr/bin/env python
#
# profile_fitting_rs_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface, Integration3DMixin

class ProfileFittingRSIntegrationExt(IntensityIface, Integration3DMixin):
  ''' Extension providing reciprocal space profile fitting. '''

  name = 'fitrs'

  phil = '''
    grid_size = 5
      .help = "The size of the reciprocal space grid for each reflection."
              "The size is the same in each dimensions"
      .type = int

    reference_frame_interval = 10
      .help = "The oscillation at which to learn new reference profiles"
      .type = int

    reference_signal_threshold = 0.02
      .help = "The threshold to use in reference profile"
      .type = float
  '''

  default = True

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithm. '''
    from dials.algorithms.integration import ProfileFittingReciprocalSpace
    self._experiments = experiments
    assert(len(experiments) == 1)
    assert(len(profile_model) == 1)
    self._algorithm = ProfileFittingReciprocalSpace(
      n_sigma = profile_model[0].n_sigma(),
      sigma_b = profile_model[0].sigma_b(deg=True),
      sigma_m = profile_model[0].sigma_m(deg=True),
      grid_size = params.integration.intensity.fitrs.grid_size,
      frame_interval = params.integration.intensity.fitrs.reference_frame_interval,
      threshold = params.integration.intensity.fitrs.reference_signal_threshold)

  def compute_intensity(self, reflections):
    ''' Compute the intensity. '''
    self._algorithm(self._experiments, reflections)
