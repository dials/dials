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

  def __init__(self, experiment, **kwargs):
    ''' Initialise a predictor for each experiment. '''
    from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
    from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor
    from dials.algorithms.spot_prediction import StillsReflectionPredictor
    from dxtbx.imageset import ImageSweep
    from dials.array_family import flex

    class Predictor(object):
      def __init__(self, name, func):
        self.name = name
        self.func = func
      def __call__(self):
        return self.func()

    # Get the force static flag
    force_static = kwargs.get("force_static", False)

    # Select the predictor class
    if isinstance(experiment.imageset, ImageSweep):
      nsp = experiment.crystal.num_scan_points
      nim = experiment.scan.get_num_images()
      if not force_static and nsp == nim + 1:
        predictor = ScanVaryingReflectionPredictor(experiment, **kwargs)
        A = [experiment.crystal.get_A_at_scan_point(i) for i in
               range(experiment.crystal.num_scan_points)]
        predict = Predictor(
          "scan varying prediction",
          lambda: predictor.for_ub(flex.mat3_double(A)))
      else:
        predictor = ScanStaticReflectionPredictor(experiment, **kwargs)
        predict = Predictor(
          "scan static prediction",
          lambda: predictor.for_ub(experiment.crystal.get_A()))
    else:
      predictor = StillsReflectionPredictor(experiment, **kwargs)

      predict = Predictor(
        "stills prediction",
        lambda: predictor.for_ub(experiment.crystal.get_A()))

    # Create and add the predictor class
    self._predict = predict

  def __call__(self):
    ''' Predict all the observable reflections.

    Returns:
      A reflection table

    '''
    from dials.util.command_line import Command
    print ' Prediction type: %s' % self._predict.name
    Command.start('Predicting reflections')
    table = self._predict()
    Command.end('Predicted %d reflections' % len(table))
    return table

  def predictor(self, index):
    ''' Get the predictor for the given experiment index. '''
    return self._predict[index]
