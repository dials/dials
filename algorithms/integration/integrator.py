#!/usr/bin/env python
#
# dials.algorithms.integration.integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Integrator(object):
  ''' The integrator base class. '''

  def __init__(self, n_sigma, n_blocks, filter_by_zeta, params):
    ''' Initialise the integrator base class.

    Params:

    '''
    self.n_sigma = n_sigma
    self.n_blocks = n_blocks
    self.filter_by_zeta = filter_by_zeta
    self.params = params

  def __call__(self, experiments, reference=None, extracted=None):
    ''' Call to integrate.

    Params:
        sweep The sweep to process
        crystal The crystal to process
        reflections The reflection list
        reference The reference profiles

    Returns:
        A reflection list

    '''
    from dials.algorithms.shoebox import ReflectionBlockExtractor
    from dials.model.data import ReflectionList
    from dials.array_family import flex
    from dials.util.command_line import Command
    from dials.algorithms.integration.lp_correction import correct_intensity

    assert(len(experiments) == 1)

    # Predict a load of reflections
    if extracted == None:
      predicted = flex.reflection_table.from_predictions(experiments)
    else:
      predicted = None

    # Get the extractor
    extract = ReflectionBlockExtractor(experiments[0], predicted,
      self.n_sigma, self.n_blocks, self.filter_by_zeta, extracted)

    # Loop through all the blocks
    result = flex.reflection_table()
    print ''
    for reflections in extract:

      reflections.integrate(experiments[0], self.params)#, reference
      result.extend(reflections)
      print ''

    # Return the reflections
    result.sort('miller_index')
    return result



class IntegratorFactory(object):
  ''' Factory class to create integrators '''

  @staticmethod
  def from_parameters(params):
    ''' Given a set of parameters, construct the integrator

    Params:
        params The input parameters

    Returns:
        The integrator instance

    '''

    # Return the integrator with the given strategies
    return Integrator(n_sigma = params.integration.shoebox.n_sigma,
                      n_blocks = params.integration.shoebox.n_blocks,
                      filter_by_zeta = params.integration.filter.by_zeta,
                      params=params)
