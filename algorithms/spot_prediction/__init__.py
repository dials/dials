from __future__ import absolute_import, division
from cctbx import sgtbx # import dependency
from dials.array_family import flex # import dependency
from dials_algorithms_spot_prediction_ext import *

# Override constructor with factory
_ScanStaticReflectionPredictor = ScanStaticReflectionPredictor

def ScanStaticReflectionPredictor(experiment, dmin=None, padding=0, **kwargs):
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
    dmin,
    padding)


# Override constructor with factory
_ScanVaryingReflectionPredictor = ScanVaryingReflectionPredictor

def ScanVaryingReflectionPredictor(experiment, dmin=None, margin=1, padding=0, **kwargs):
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
    margin,
    padding)


def StillsReflectionPredictor(experiment, dmin=None, spherical_relp=False,
                              **kwargs):
  '''
  A factory function for the reflection predictor.

  :param experiment: The experiment to predict for
  :param dmin: The maximum resolution to predict to
  :param spherical_relp: Whether to use the spherical relp prediction model
  :return: The spot predictor

  '''

  # FIXME Selection of reflection predictor type is ugly. What is a better
  # way here? Should it be based entirely on the existence of certain types
  # of profile model within the experiment?

  # Get dmin if it is not set
  if dmin is None:
    dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

  if spherical_relp:
    return SphericalRelpStillsReflectionPredictor(
      experiment.beam,
      experiment.detector,
      experiment.crystal.get_A(),
      experiment.crystal.get_unit_cell(),
      experiment.crystal.get_space_group().type(),
      dmin)

  # Create the reflection predictor
  try:
    if experiment.crystal.get_half_mosaicity_deg() is not None and experiment.crystal.get_domain_size_ang() is not None:
      return NaveStillsReflectionPredictor(
        experiment.beam,
        experiment.detector,
        experiment.crystal.get_A(),
        experiment.crystal.get_unit_cell(),
        experiment.crystal.get_space_group().type(),
        dmin,
        experiment.crystal.get_half_mosaicity_deg(),
        experiment.crystal.get_domain_size_ang())
  except AttributeError:
    pass

  return StillsDeltaPsiReflectionPredictor(
    experiment.beam,
    experiment.detector,
    experiment.crystal.get_A(),
    experiment.crystal.get_unit_cell(),
    experiment.crystal.get_space_group().type(),
    dmin)
