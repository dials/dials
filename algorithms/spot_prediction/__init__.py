from __future__ import division
from cctbx import sgtbx # import dependency
from dials.array_family import flex # import dependency
from dials.model.data import Reflection, ReflectionList # import dependency
from dials_algorithms_spot_prediction_ext import *

_ScanStaticReflectionPredictor = ScanStaticReflectionPredictor

def ScanStaticReflectionPredictor(experiment, dmin=None):
  ''' A constructor for the reflection predictor. '''

  # Get dmin if it is not set
  if dmin is None:
    dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

  # Create the reflection predictor
  return _ScanStaticReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.goniometer,
    experiment.scan,
    experiment.crystal.get_unit_cell(),
    experiment.crystal.get_space_group().type(),
    experiment.crystal.get_A(),
    dmin)
