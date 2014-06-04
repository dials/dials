#
# integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Integrator(object):
  ''' A class to perform the integration. '''

  def __init__(self, params, experiments, extractor=None):
    ''' Initialise the integrator. '''
    from math import pi
    from dials.algorithms import shoebox

    # Ensure we have 1 experiment at the moment
    assert(len(experiments) == 1)
    assert(extractor is not None)

    # Save the parameters
    self.params = params
    self.experiments = experiments
    self.extractor = extractor

    # Create the shoebox masker
    n_sigma = params.integration.shoebox.n_sigma
    sigma_b = params.integration.shoebox.sigma_b * pi / 180.0
    sigma_m = params.integration.shoebox.sigma_m * pi / 180.0
    delta_b = n_sigma * sigma_b
    delta_m = n_sigma * sigma_m
    self._mask_profiles = shoebox.Masker(experiments[0], delta_b, delta_m)

  def integrate(self):
    ''' Integrate all the reflections. '''
    from dials.array_family import flex
    result = flex.reflection_table()
    for indices, reflections in self.extractor:
      self._mask_profiles(reflections, None)
      reflections.integrate(self.experiments[0])
      del reflections['shoebox']
      del reflections['rs_shoebox']
      result.extend(reflections)
    assert(len(result) > 0)
    result.sort('miller_index')
    return result
