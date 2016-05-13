#!/usr/bin/env python
#
# simple_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class SimpleBackgroundExt(BackgroundIface):
  ''' An extension class implementing simple background subtraction. '''

  name = 'simple'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''
      outlier
        .help = "Outlier rejection prior to background fit"
      {
        algorithm = *null nsigma truncated normal mosflm tukey
          .help = "The outlier rejection algorithm."
          .type = choice

        nsigma
          .help = "Parameters for nsigma outlier rejector"
          .expert_level = 1
        {
          lower = 3
            .help = "Lower n sigma"
            .type = float
          upper = 3
            .help = "Upper n sigma"
            .type = float
        }

        truncated
          .help = "Parameters for truncated outlier rejector"
          .expert_level = 1
        {
          lower = 0.01
            .help = "Lower bound"
            .type = float
          upper = 0.01
            .help = "Upper bound"
            .type = float
        }

        normal
          .help = "Parameters for normal outlier rejector"
          .expert_level = 1
        {
          min_pixels = 10
            .help = "The minimum number of pixels to use in calculating the"
                    "background intensity."
            .type = int
        }

        mosflm
          .help = "Parameters for mosflm-like outlier rejector. This algorithm"
                  "is mainly used in conjunction with a linear 2d background."
          .expert_level = 1
        {
          fraction = 1.0
            .help = "The fraction of pixels to use in determining the initial"
                    "plane used for outlier rejection."
            .type = float

          n_sigma = 4.0
            .help = "The number of standard deviations above the threshold plane"
                    "to use in rejecting outliers from background calculation."
            .type = float
        }

        tukey
          .help = "Parameters for tukey outlier rejector"
          .expert_level = 1
        {
          lower = 1.5
            .help = "Lower IQR multiplier"
            .type = float
          upper = 1.5
            .help = "Upper IQR multiplier"
            .type = float
        }
      }

      model
        .help = "Background model"
      {
        algorithm = constant2d *constant3d linear2d linear3d
          .help = "The choice of background model"
          .type = choice
      }

    ''')
    return phil

  def __init__(self, params, experiments):
    '''
    Initialise the algorithm.

    :param params: The input parameters
    :param experiments: The list of experiments

    '''
    from libtbx.phil import parse
    from dials.algorithms.background.simple import BackgroundAlgorithm

    # Create some default parameters
    if params is None:
      params = self.phil().fetch(parse('')).extract()
    else:
      params = params.integration.background.simple

    # Create some keyword parameters
    kwargs = {
      'model' : params.model.algorithm,
      'outlier' : params.outlier.algorithm
    }

    # Create all the keyword parameters
    if params.outlier.algorithm == 'null':
      pass
    elif params.outlier.algorithm == 'truncated':
      kwargs['lower'] = params.outlier.truncated.lower
      kwargs['upper'] = params.outlier.truncated.upper
    elif params.outlier.algorithm == 'nsigma':
      kwargs['lower'] = params.outlier.nsigma.lower
      kwargs['upper'] = params.outlier.nsigma.upper
    elif params.outlier.algorithm == 'normal':
      kwargs['min_pixels'] = params.outlier.normal.min_pixels
    elif params.outlier.algorithm == 'mosflm':
      kwargs['fraction'] = params.outlier.mosflm.fraction
      kwargs['n_sigma'] = params.outlier.mosflm.n_sigma
    elif params.outlier.algorithm == 'tukey':
      kwargs['lower'] = params.outlier.tukey.lower
      kwargs['upper'] = params.outlier.tukey.upper

    # Create the algorithm
    self._algorithm = BackgroundAlgorithm(experiments, **kwargs)

  def compute_background(self, reflections, image_volume=None):
    '''
    Compute the background.

    :param reflections: The list of reflections

    '''
    return self._algorithm.compute_background(
      reflections, image_volume=image_volume)
