from __future__ import division
from dials.array_family import flex # import dependency
from dials_algorithms_profile_model_gaussian_rs_ext import *

from model import phil_scope        # implicit dependency
from model import ProfileModel      # implicit dependency
from model import ProfileModelList  # implicit dependency


def BBoxCalculator(experiment, delta_b, delta_m):
  ''' Return the relavent bbox calculator. '''
  if experiment.goniometer is None or experiment.scan is None:
    algorithm = BBoxCalculator2D(
      experiment.beam,
      experiment.detector,
      delta_b,
      delta_m)
  else:
    algorithm = BBoxCalculator3D(
      experiment.beam,
      experiment.detector,
      experiment.goniometer,
      experiment.scan,
      delta_b,
      delta_m)
  return algorithm
