from __future__ import division
from cctbx import sgtbx # import dependency
from dials.array_family import flex # import dependency
from dials.model.data import Reflection, ReflectionList # import dependency
from dials_algorithms_spot_prediction_ext import *

# Override constructor with factory
_ScanStaticReflectionPredictor = ScanStaticReflectionPredictor

def ScanStaticReflectionPredictor(experiment, dmin=None, **kwargs):
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
    dmin)


# Override constructor with factory
_ScanVaryingReflectionPredictor = ScanVaryingReflectionPredictor

def ScanVaryingReflectionPredictor(experiment, dmin=None, margin=1, **kwargs):
  ''' A constructor for the reflection predictor. '''
  from dials.array_family import flex

  # Get dmin if it is not set
  if dmin is None:
    dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

  # Create the reflection predictor
  return _ScanVaryingReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.goniometer,
    experiment.scan,
    dmin,
    margin)


# Override constructor with factory
_StillsReflectionPredictor = StillsReflectionPredictor

def StillsReflectionPredictor(experiment, **kwargs):
  ''' A constructor for the reflection predictor. '''

  # Create the reflection predictor
  return _StillsReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.crystal.get_A())
