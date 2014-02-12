#!/usr/bin/env python
#
# reflection_predictor.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ReflectionPredictor(object):
  ''' A reflection predictor that takes a number of experiments and does the
  proper prediction for each type of experiment. '''

  def __init__(self, experiments, **kwargs):
    ''' Initialise a predictor for each experiment. '''
    from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
    from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor
    from dials.algorithms.spot_prediction import StillsReflectionPredictor
    from dxtbx.imageset import ImageSweep

    # Create all the reflection predictors
    self._predict = []
    for e in experiments:
      # Select the predictor class
      if isinstance(e.imageset, ImageSweep):
        if e.crystal.num_scan_points == e.scan.get_num_images() + 1:
          Predictor = ScanVaryingReflectionPredictor
        else:
          Predictor = ScanStaticReflectionPredictor
      else:
        Predictor = StillsReflectionPredictor

      # Create and add the predictor class
      self._predict.append(Predictor(e, **kwargs))

  def __call__(self):
    ''' Predict all the observable reflections.

    Returns:
      A reflection table

    '''
    from dials.array_family import flex
    from dials.util.command_line import Command
    table = flex.reflection_table()
    for i, predict in enumerate(self._predict):
      Command.start('Predicting reflections for experiment %d' % i)
      temp = predict()
      temp['id'] = flex.size_t(temp.nrows(), i)
      table.extend(temp)
      Command.end('Predicted %d reflections for experiment %s' % (len(temp), i))
    return table

  def predictor(self, index):
    ''' Get the predictor for the given experiment index. '''
    return self._predict[index]

