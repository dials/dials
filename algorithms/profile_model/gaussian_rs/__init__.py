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


def PartialityCalculator(experiment, delta_m):
  ''' Return the relavent partiality calculator. '''
  if experiment.goniometer is None or experiment.scan is None:
    print "WARNING: Stills partiality is currently a placeholder"
    algorithm = PartialityCalculator2D(
      experiment.beam,
      delta_m)
  else:
    algorithm = PartialityCalculator3D(
      experiment.beam,
      experiment.goniometer,
      experiment.scan,
      delta_m)
  return algorithm


def MaskCalculator(experiment, delta_b, delta_m):
  ''' Return the relavent partiality calculator. '''
  if experiment.goniometer is None or experiment.scan is None:
    algorithm = MaskCalculator2D(
      experiment.beam,
      experiment.detector,
      delta_b,
      delta_m)
  else:
    algorithm = MaskCalculator3D(
      experiment.beam,
      experiment.detector,
      experiment.goniometer,
      experiment.scan,
      delta_b,
      delta_m)
  return algorithm
