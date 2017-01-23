#!/usr/bin/env python
#
# general.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

class BackgroundAlgorithm(object):
  ''' Class to do background subtraction. '''

  def __init__(self, experiments,
               outlier='nsigma',
               model='constant3d',
               **kwargs):
    '''
    Initialise the algorithm.

    :param experiments: The list of experiments
    :param outlier: The outlier rejection algorithm
    :param model: The background model algorithm

    '''
    from dials.algorithms.background.simple import Creator
    from dials.algorithms.background.simple import TruncatedOutlierRejector
    from dials.algorithms.background.simple import NSigmaOutlierRejector
    from dials.algorithms.background.simple import NormalOutlierRejector
    from dials.algorithms.background.simple import MosflmOutlierRejector
    from dials.algorithms.background.simple import TukeyOutlierRejector
    from dials.algorithms.background.simple import Constant2dModeller
    from dials.algorithms.background.simple import Constant3dModeller
    from dials.algorithms.background.simple import Linear2dModeller
    from dials.algorithms.background.simple import Linear3dModeller

    def select_modeller():
      if model == 'constant2d':
        return Constant2dModeller()
      elif model == 'constant3d':
        return Constant3dModeller()
      elif model == 'linear2d':
        return Linear2dModeller()
      elif model == 'linear3d':
        return Linear3dModeller()
      raise RuntimeError("Unexpected background model: %s" % model)

    def select_rejector():
      if outlier == 'null':
        return None
      elif outlier == 'truncated':
        return TruncatedOutlierRejector(
          kwargs.get("lower", 0.01),
          kwargs.get("upper", 0.01))
      elif outlier == 'nsigma':
        return NSigmaOutlierRejector(
          kwargs.get("lower", 3),
          kwargs.get("upper", 3))
      elif outlier == 'normal':
        return NormalOutlierRejector(
          kwargs.get("min_pixels", 10))
      elif outlier == 'plane':
        return MosflmOutlierRejector(
          kwargs.get("fraction", 1.0),
          kwargs.get("n_sigma", 4.0))
      elif outlier == 'tukey':
        return TukeyOutlierRejector(
          kwargs.get("lower", 1.5),
          kwargs.get("upper", 1.5))
      raise RuntimeError("Unexpected outlier rejector: %s" % outlier)

    # Get the minimum number of pixels
    min_pixels = kwargs.get("min_pixels", 10)

    modeller = select_modeller()
    rejector = select_rejector()
    self._creator = Creator(modeller, rejector, min_pixels=min_pixels)

  def compute_background(self, reflections, image_volume=None):
    '''
    Compute the backgrond.

    :param reflections: The list of reflections

    '''
    from dials.array_family import flex

    # Do the background subtraction
    if image_volume is None:
      reflections['background.mse'] = flex.double(len(reflections))
      reflections['background.dispersion'] = flex.double(len(reflections))
      success = self._creator(
        reflections['shoebox'],
        reflections['background.mse'],
        reflections['background.dispersion'])
      reflections['background.mean'] = reflections['shoebox'].mean_background()
    else:
      success = self._creator(reflections, image_volume)
    reflections.set_flags(success != True, reflections.flags.dont_integrate)
    return success
