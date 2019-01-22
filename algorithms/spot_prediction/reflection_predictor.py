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

from __future__ import absolute_import, division
import logging
logger = logging.getLogger(__name__)

from libtbx.phil import parse
from libtbx.utils import Sorry

# The phil parameters
phil_scope = parse('''

  prediction {

    d_min = None
      .type = float
      .help = "The maximum resolution limit"

    d_max = None
      .type = float
      .help = "The minimum resolution limit"

    margin = 1
      .type = int
      .help = "The margin to use to scan varying prediction"

    force_static = False
      .type = bool
      .help = "For scan-varying prediction force scan-static prediction"

    padding = 1.0
      .type = float(value_min=0)
      .help = "The padding in degrees"
  }

''')


class ReflectionPredictor(object):
  '''
  A reflection predictor that takes a number of experiments and does the proper
  prediction for each type of experiment.

  '''

  def __init__(self,
               experiment,
               dmin=None,
               dmax=None,
               margin=1,
               force_static=False,
               padding=0):
    '''
    Initialise a predictor for each experiment.

    :param experiment: The experiment to predict for
    :param dmin: The maximum resolution
    :param dmax: The minimum resolution
    :param margin: The margin of hkl to predict
    :param force_static: force scan varying prediction to be static

    '''
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
        result = self.func()
        if dmax is not None:
          assert(dmax > 0)
          result.compute_d_single(experiment)
          mask = result['d'] > dmax
          result.del_selected(mask)
        return result

    # Check prediction to maximum resolution is possible
    wl = experiment.beam.get_wavelength()
    if dmin is not None and dmin < 0.5 * wl:
      raise Sorry("Prediction at d_min of {0} is not possible "
                  "with wavelength {1}".format(dmin, wl))

    # Select the predictor class
    if isinstance(experiment.imageset, ImageSweep):
      xl_nsp = experiment.crystal.num_scan_points
      bm_nsp = experiment.beam.num_scan_points
      gn_nsp = experiment.goniometer.num_scan_points
      nim = experiment.scan.get_num_images()

      sv_compatible = (xl_nsp == nim + 1) or (bm_nsp == nim + 1)
      if not force_static and sv_compatible:
        predictor = ScanVaryingReflectionPredictor(
          experiment,
          dmin=dmin,
          margin=margin,
          padding=padding)

        if bm_nsp == 0 and gn_nsp == 0:
          # Only varying crystal
          A = [experiment.crystal.get_A_at_scan_point(i) for i in
                 range(experiment.crystal.num_scan_points)]
          predict = Predictor(
            "scan varying crystal prediction",
            lambda: predictor.for_ub(flex.mat3_double(A)))

        else:
          # Any allowed model may vary
          if xl_nsp == nim + 1:
            A = [experiment.crystal.get_A_at_scan_point(i) for i in
                 range(experiment.crystal.num_scan_points)]
          else:
            A = [experiment.crystal.get_A() for i in range(nim + 1)]
          if bm_nsp == nim + 1:
            s0 = [experiment.beam.get_s0_at_scan_point(i) for i in
                 range(experiment.beam.num_scan_points)]
          else:
            s0 = [experiment.beam.get_s0() for i in range(nim + 1)]
          if gn_nsp == nim + 1:
            S =  [experiment.goniometer.get_setting_rotation_at_scan_point(i)
                  for i in range(experiment.goniometer.num_scan_points)]
          else:
            S = [experiment.goniometer.get_setting_rotation()
                 for i in range(nim + 1)]
          predict = Predictor(
            "scan varying model prediction",
            lambda: predictor.for_varying_models(
                flex.mat3_double(A),
                flex.vec3_double(s0),
                flex.mat3_double(S)))
      else:
        predictor = ScanStaticReflectionPredictor(
          experiment,
          dmin=dmin,
          padding=padding)

        # Choose index generation method based on number of images
        # https://github.com/dials/dials/issues/585
        if experiment.scan.get_num_images() > 50:
          predict_method = predictor.for_ub_old_index_generator
        else:
          predict_method = predictor.for_ub

        predict = Predictor(
          "scan static prediction",
          lambda: predict_method(experiment.crystal.get_A()))
    else:
      predictor = StillsReflectionPredictor(
        experiment,
        dmin=dmin)

      predict = Predictor(
        "stills prediction",
        lambda: predictor.for_ub(experiment.crystal.get_A()))

    # Create and add the predictor class
    self._predict = predict

  def __call__(self):
    '''
    Predict all the observable reflections.

    :return: A reflection table

    '''
    logger.info('Prediction type: %s' % self._predict.name)
    table = self._predict()
    logger.info('Predicted %d reflections' % len(table))
    return table

  def predictor(self, index):
    '''
    Get the predictor for the given experiment index.

    :param index: The experiment index
    :return: The predictor

    '''
    return self._predict[index]
