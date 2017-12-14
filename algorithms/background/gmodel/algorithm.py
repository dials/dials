#!/usr/bin/env python
#
# algorithm.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

class ModelCache(object):
  '''
  A class to cache the model

  '''
  def __init__(self):
    '''
    Create a model dictionary

    '''
    self.model = dict()

  def get(self, name):
    '''
    Get the model

    '''
    if name is None:
      raise RuntimeError('Model is not specified')
    try:
      model = self.model[name]
    except KeyError:
      import cPickle as pickle
      with open(name) as infile:
        model = pickle.load(infile)
        self.model[name] = model
    return model


# Instance of the model cache
global_model_cache = ModelCache()


class BackgroundAlgorithm(object):
  ''' Class to do background subtraction. '''

  def __init__(self,
               experiments,
               model=None,
               robust=False,
               tuning_constant=1.345,
               min_pixels=10):
    '''
    Initialise the algorithm.

    :param experiments: The list of experiments
    :param model: The background model
    :param robust: Use the robust background algorithm
    :param tuning_constant: The robust tuning constant

    '''
    from dials.algorithms.background.gmodel import Creator

    # Get the model
    model = global_model_cache.get(model)

    # Create the background creator
    self._create = Creator(
      model           = model,
      robust          = robust,
      tuning_constant = tuning_constant,
      max_iter        = 100,
      min_pixels      = min_pixels)

  def compute_background(self, reflections, image_volume=None):
    '''
    Compute the backgrond.

    :param reflections: The list of reflections

    '''
    from dials.array_family import flex

    # Do the background subtraction
    if image_volume is None:
      success = self._create(reflections)
      reflections['background.mean'] = reflections['shoebox'].mean_modelled_background()
    else:
      success = self._create(reflections, image_volume)
    reflections.set_flags(success != True, reflections.flags.dont_integrate)
    return success


class GModelBackgroundCalculatorFactory(object):
  ''' Class to do background subtraction. '''

  @classmethod
  def create(Class,
             experiments,
             model=None,
             robust=False,
             tuning_constant=1.345,
             min_pixels=10):
    '''
    Initialise the algorithm.

    :param experiments: The list of experiments
    :param model: The background model
    :param robust: Use the robust background algorithm
    :param tuning_constant: The robust tuning constant

    '''
    from dials.algorithms.integration.parallel_integrator import GModelBackgroundCalculator
    from dials.algorithms.background.gmodel import Creator

    # Get the model
    model = global_model_cache.get(model)

    # Create the background creator
    return GModelBackgroundCalculator(
      model           = model,
      robust          = robust,
      tuning_constant = tuning_constant,
      max_iter        = 100,
      min_pixels      = min_pixels)
