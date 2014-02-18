#!/usr/bin/env python
#
# dials.algorithms.integration.background_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class BackgroundFactory(object):
  ''' Factory class to create integrators '''

  @staticmethod
  def from_parameters(params):
    ''' Given a set of parameters, construct the integrator

    Params:
        params The input parameters

    Returns:
        The background calculator instance

    '''
    from dials.algorithms.background import NullSubtractor
    from dials.algorithms.background import XdsSubtractor
    from dials.algorithms.background import FableSubtractor
    from dials.algorithms.background import FlatSubtractor
    from dials.algorithms.background import CurvedSubtractor
    from dials.algorithms.background import InclinedSubtractor

    # Shorten parameter path
    integration = params.integration

    # Configure the NULL subtractor
    if (integration.background.algorithm == 'none' or
        integration.background.algorithm == None):
      algorithm = NullSubtractor()

    # Configure the XDS subtractor
    elif integration.background.algorithm == 'xds':
      algorithm = XdsSubtractor(
          min_data = integration.background.min_pixels,
          n_sigma = integration.background.n_sigma)

    # Configure the Fable subtractor
    elif params.integration.background.algorithm == 'fable':
      algorithm = FableSubtractor(
          min_data = integration.background.min_pixels,
          n_sigma = integration.background.n_sigma)

    # Configure the flat subtractor
    elif integration.background.algorithm == 'flat':
      algorithm = FlatSubtractor()

    # Configure the inclined plane subtractor
    elif integration.background.algorithm == 'inclined':
      algorithm = InclinedSubtractor()

    # Configure the esmerelda curved subtractor
    elif integration.background.algorithm == 'esmeralda':
      algorithm = CurvedSubtractor()

    # Unknown subtractor
    else:
      raise RuntimeError('Unknown background algorithm')

    # Return the algorithm
    return algorithm
