#
# parallel_integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
from dials_algorithms_integration_parallel_integrator_ext import *

from dials.algorithms.background.simple.algorithm import SimpleBackgroundCalculatorFactory
from dials.algorithms.background.glm.algorithm import GLMBackgroundCalculatorFactory
from dials.algorithms.background.gmodel.algorithm import GModelBackgroundCalculatorFactory


class GaussianRSIntensityCalculatorFactory(object):
  '''
  A class to create the intensity calculator

  '''

  @classmethod
  def create(Class,
             data,
             detector_space = False,
             deconvolution = False):
    '''
    Create the intensity calculator

    '''

    # Return the intensity algorithm
    return GaussianRSIntensityCalculator(data, detector_space, deconvolution)
