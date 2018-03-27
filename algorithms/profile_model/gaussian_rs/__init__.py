from __future__ import absolute_import, division
from __future__ import print_function
from dials.array_family import flex # import dependency
from dials.algorithms.profile_model import modeller # import dependency
from dials_algorithms_profile_model_gaussian_rs_ext import *

from dials.algorithms.profile_model.gaussian_rs.model import phil_scope        # implicit dependency
from dials.algorithms.profile_model.gaussian_rs.model import Model             # implicit dependency


def BBoxCalculator(crystal, beam, detector, goniometer, scan, delta_b, delta_m):
  ''' Return the relavent bbox calculator. '''
  if goniometer is None or scan is None or scan.get_oscillation()[1] == 0:
    algorithm = BBoxCalculator2D(
      beam,
      detector,
      delta_b,
      delta_m)
  else:
    algorithm = BBoxCalculator3D(
      beam,
      detector,
      goniometer,
      scan,
      delta_b,
      delta_m)
  return algorithm


def PartialityCalculator(crystal, beam, detector, goniometer, scan, sigma_m):
  ''' Return the relavent partiality calculator. '''
  if goniometer is None or scan is None or scan.get_oscillation()[1] == 0:
    print("WARNING: Stills partiality is currently a placeholder")
    algorithm = PartialityCalculator2D(
      beam,
      sigma_m)
  else:
    algorithm = PartialityCalculator3D(
      beam,
      goniometer,
      scan,
      sigma_m)
  return algorithm


def MaskCalculator(crystal, beam, detector, goniometer, scan, delta_b, delta_m):
  ''' Return the relavent partiality calculator. '''
  if goniometer is None or scan is None or scan.get_oscillation()[1] == 0:
    algorithm = MaskCalculator2D(
      beam,
      detector,
      delta_b,
      delta_m)
  else:
    algorithm = MaskCalculator3D(
      beam,
      detector,
      goniometer,
      scan,
      delta_b,
      delta_m)
  return algorithm
