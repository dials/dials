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

  def __init__(self, experiments):
    ''' Initialise a predictor for each experiment. '''
    from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
    from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor
    from dials.algorithms.spot_prediction import ScanStillsReflectionPredictor

    # Create all the reflection predictors
    self._predict = []
    for e in experiments:

      # Select the predictor class
      if isinstance(e.imageset, ImageSweep):
        if e.crystal.is_scan_varying():
          Predictor = ScanVaryingReflectionPredictor
        else:
          Predictor = ScanStaticReflectionPredictor
      else:
        Predictor = StillsReflectionPredictor

      # Create and add the predictor class
      self._predict.append(Predictor(beam, detector,
        goniometer, scan, crystal))

  def all_possible(self):
    ''' Predict all the observable reflections.

    Returns:
      A reflection table

    '''
    from dials.array_family import flex
    table = flex.reflection_table()
    for i, predict in enumerate(self._predict):
      temp = predict.all_possible()
      temp['id'] = flex.size_t(temp.nrows(), i)
      table.extend(temp)
    return table

  def selected(self, hkl, experiment, panel=None):
    ''' Predict the given reflections.

    Params:
      hkl The miller indices to predict for.
      experiment The experiment of each hkl.
      panel The panel number of each hkl (optional)

    Returns
      A reflection table

    '''
    from dials.array_family import flex
    table = flex.reflection_table()
    for i, predict in enumerate(self._predict):
      mask = experiment == i
      if panel:
        temp = predict.selected(hkl.select(mask), panel.select(mask))
      else:
        temp = predict.selected(hkl.select(mask))
      temp['id'] = flex.size_t(temp.nrows(), i)
      table.extend(temp)
    return table
