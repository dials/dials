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


class MaskCalculatorFactory(object):
  '''
  A factory function to return a mask calculator object

  '''

  @classmethod
  def create(Class, experiments, params=None):
    '''
    Select the mask calculator

    '''
    from dials.algorithms.profile_model.gaussian_rs.algorithm import GaussianRSMaskCalculatorFactory

    # Get the parameters
    if params is None:
      from dials.command_line.integrate import phil_scope
      params = phil_scope.extract()

    # Select the factory function
    selection = params.profile.algorithm
    if selection == "gaussian_rs":
      algorithm = GaussianRSMaskCalculatorFactory.create(experiments)
    else:
      raise RuntimeError("Unknown profile model algorithm")

    # Create the mask algorithm
    return algorithm


class BackgroundCalculatorFactory(object):
  '''
  A factory function to return a background calculator object

  '''

  @classmethod
  def create(Class, experiments, params=None):
    '''
    Select the background calculator

    '''
    from dials.algorithms.background.simple.algorithm import SimpleBackgroundCalculatorFactory
    from dials.algorithms.background.glm.algorithm import GLMBackgroundCalculatorFactory
    from dials.algorithms.background.gmodel.algorithm import GModelBackgroundCalculatorFactory

    # Get the parameters
    if params is None:
      from dials.command_line.integrate import phil_scope
      params = phil_scope.extract()

    # Select the factory function
    selection = params.integration.background.algorithm
    if selection == "simple":

      # Get parameters
      params = params.integration.background.simple

      # Create some keyword parameters
      kwargs = {
        'model'      : params.model.algorithm,
        'outlier'    : params.outlier.algorithm,
        'min_pixels' : params.min_pixels
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
      elif params.outlier.algorithm == 'plane':
        kwargs['fraction'] = params.outlier.plane.fraction
        kwargs['n_sigma'] = params.outlier.plane.n_sigma
      elif params.outlier.algorithm == 'tukey':
        kwargs['lower'] = params.outlier.tukey.lower
        kwargs['upper'] = params.outlier.tukey.upper

      # Create the algorithm
      algorithm = SimpleBackgroundCalculatorFactory.create(experiments, **kwargs)

    elif selection == "glm":

      # Get the parameters
      params = params.integration.background.glm

      # Create the algorithm
      algorithm = GLMBackgroundCalculatorFactory.create(
        experiments,
        model           = params.model.algorithm,
        tuning_constant = params.robust.tuning_constant,
        min_pixels      = params.min_pixels)

    elif selection == "gmodel":

      # Get the parameters
      params = params.integration.background.gmodel

      # Create the algorithm
      algorithm = GModelBackgroundCalculatorFactory.create(
        experiments,
        model           = params.model,
        robust          = params.robust.algorithm,
        tuning_constant = params.robust.tuning_constant,
        min_pixels      = params.min_pixels)

    else:
      raise RuntimeError("Unknown background algorithm")

    # Return the background calculator
    return algorithm


class IntensityCalculatorFactory(object):
  '''
  A factory function to return an intensity calculator object

  '''

  @classmethod
  def create(Class, experiments, reference_profiles, params=None):
    '''
    Select the mask calculator

    '''

    # Get the parameters
    if params is None:
      from dials.command_line.integrate import phil_scope
      params = phil_scope.extract()

    # Select the factory function
    selection = params.profile.algorithm
    if selection == "gaussian_rs":

      # Get the parameters
      params = params.profile.gaussian_rs.fitting

      # Check for detector space
      if params.fit_method == "reciprocal_space":
        detector_space = False
      elif params.fit_method == "detector_space":
        detector_space = True
      else:
        raise RuntimeError("Unknown fit method: %s" % params.fit_method)

      # Create the algorithm
      algorithm = GaussianRSIntensityCalculatorFactory.create(
        reference_profiles,
        detector_space = detector_space,
        deconvolution  = params.detector_space.deconvolution)

    else:
      raise RuntimeError("Unknown profile model algorithm")

    # Return the algorithm
    return algorithm
