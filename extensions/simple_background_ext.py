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
  ''' An extension class implementing XDS background subtraction. '''

  name = 'simple'

  default = True

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''
      outlier
        .help = "Outlier rejection prior to background fit"
      {
        algorithm = null *nsigma truncated normal mosflm
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
    ''' Initialise the algorithm. '''
    from libtbx.phil import parse
    from dials.algorithms.background.simple import BackgroundAlgorithm

    # Create some default parameters
    if params is None:
      phil = '''
        integration {
          background {
            simple {
              %s
            }
          }
        }''' % self.phil
      params = parse(phil).extract()

    # Get the model and outlier algorithms
    model = params.integration.background.simple.model
    outlier = params.integration.background.simple.outlier

    # Create some keyword parameters
    kwargs = {
      'model' : model.algorithm,
      'outlier' : outlier.algorithm
    }

    # Create all the keyword parameters
    if outlier.algorithm == 'null':
      pass
    elif outlier.algorithm == 'truncated':
      kwargs['lower'] = outlier.truncated.lower
      kwargs['upper'] = outlier.truncated.upper
    elif outlier.algorithm == 'nsigma':
      kwargs['lower'] = outlier.nsigma.lower
      kwargs['upper'] = outlier.nsigma.upper
    elif outlier.algorithm == 'normal':
      kwargs['min_pixels'] = outlier.normal.min_pixels
    elif outlier.algorithm == 'mosflm':
      kwargs['fraction'] = outlier.mosflm.fraction
      kwargs['n_sigma'] = outlier.mosflm.n_sigma

    # Create the algorithm
    self._algorithm = BackgroundAlgorithm(experiments, **kwargs)

  def compute_background(self, reflections):
    ''' Compute the backgrond. '''
    self._algorithm.compute_background(reflections)
