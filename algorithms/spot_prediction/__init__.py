from __future__ import division
from cctbx import sgtbx # import dependency
from dials.array_family import flex # import dependency
from dials_algorithms_spot_prediction_ext import *

# Override constructor with factory
_ScanStaticReflectionPredictor = ScanStaticReflectionPredictor

def ScanStaticReflectionPredictor(experiment, dmin=None, **kwargs):
  '''
  A constructor for the reflection predictor.

  :param experiment: The experiment to predict for
  :param dmin: The maximum resolution to predict to
  :return: The spot predictor

  '''

  # Get dmin if it is not set
  if dmin is None:
    dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

  # Only remove certain systematic absences
  space_group = experiment.crystal.get_space_group()
  space_group = space_group.build_derived_patterson_group()

  # Create the reflection predictor
  return _ScanStaticReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.goniometer,
    experiment.scan,
    experiment.crystal.get_unit_cell(),
    space_group.type(),
    dmin)


# Override constructor with factory
_ScanVaryingReflectionPredictor = ScanVaryingReflectionPredictor

def ScanVaryingReflectionPredictor(experiment, dmin=None, margin=1, **kwargs):
  '''
  A constructor for the reflection predictor.

  :param experiment: The experiment to predict for
  :param dmin: The maximum resolution to predict to
  :param margin: The margin for prediction
  :return: The spot predictor

  '''

  # Get dmin if it is not set
  if dmin is None:
    dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

  # Only remove certain systematic absences
  space_group = experiment.crystal.get_space_group()
  space_group = space_group.build_derived_patterson_group()

  # Create the reflection predictor
  return _ScanVaryingReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.goniometer,
    experiment.scan,
    space_group.type(),
    dmin,
    margin)


# Override constructor with factory
_StillsReflectionPredictor = StillsDeltaPsiReflectionPredictor

def StillsReflectionPredictor(experiment, dmin=None, **kwargs):
  '''
  A constructor for the reflection predictor.

  :param experiment: The experiment to predict for
  :param dmin: The maximum resolution to predict to
  :return: The spot predictor

  '''

  # Get dmin if it is not set
  if dmin is None:
    dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

  # Create the reflection predictor
  return _StillsReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.crystal.get_A(),
    experiment.crystal.get_unit_cell(),
    experiment.crystal.get_space_group().type(),
    dmin)
