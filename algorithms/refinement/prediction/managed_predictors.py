#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Managed reflection prediction for refinement.

* ScansRayPredictor adapts DIALS prediction for use in refinement, by keeping
  up to date with the current model geometry
* StillsRayPredictor predicts reflections without a goniometer, under
  the naive assumption that the relp is already in reflecting position

"""

from __future__ import absolute_import, division

from math import pi
from dials.algorithms.spot_prediction import ScanStaticRayPredictor

from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor as sc
from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor as sv
from dials.algorithms.spot_prediction import StillsReflectionPredictor as st

class ScansRayPredictor(object):
  """
  Predict for a relp based on the current states of models of the
  experimental geometry. This is a wrapper for DIALS' C++
  RayPredictor class, which does the real work. This class keeps track
  of the experimental geometry, and instantiates a RayPredictor when
  required.
  """

  def __init__(self, experiments, sweep_range=(0, 2.*pi)):
    """Construct by linking to instances of experimental model classes"""

    self._experiments = experiments
    self._sweep_range = sweep_range

  def __call__(self, hkl, experiment_id=0, UB=None):
    """
    Solve the prediction formula for the reflecting angle phi.

    If UB is given, override the contained crystal model. This is
    for use in refinement with time-varying crystal parameters
    """

    e = self._experiments[experiment_id]
    ray_predictor = ScanStaticRayPredictor(
      e.beam.get_s0(),
      e.goniometer.get_rotation_axis_datum(),
      e.goniometer.get_fixed_rotation(),
      e.goniometer.get_setting_rotation(),
      self._sweep_range)

    UB_ = UB if UB else e.crystal.get_A()

    rays = ray_predictor(hkl, UB_)

    return rays

class ExperimentsPredictor(object):
  """
  Predict for relps based on the current states of models of the experimental
  geometry. This version manages multiple experiments, selecting the correct
  predictor in each case.
  """

  def __init__(self, experiments, force_stills=False, spherical_relp=False):
    """Construct by linking to instances of experimental model classes"""

    self._experiments = experiments
    self._force_stills = force_stills
    self._spherical_relp = spherical_relp

  def __call__(self, reflections):
    """Predict for all reflections at the current model geometry"""

    if self._force_stills:
      predictors = [st(e, spherical_relp=self._spherical_relp) \
                          for e in self._experiments]
    else:
      predictors = [sc(e) if e.goniometer else st(e,
        spherical_relp=self._spherical_relp) for e in self._experiments]
    self._UBs = [e.crystal.get_U() * e.crystal.get_B() for e in self._experiments]

    for iexp, e in enumerate(self._experiments):

      # select the reflections for this experiment only
      sel = reflections['id'] == iexp
      refs = reflections.select(sel)

      # stills
      if not e.goniometer or self._force_stills:
        predictor = st(e, spherical_relp=self._spherical_relp)
        UB = e.crystal.get_A()
        predictor.for_reflection_table(refs, UB)
      # scan-varying
      elif 'ub_matrix' in refs:
        predictor = sv(e)
        UB = refs['ub_matrix']
        s0 = refs['s0_vector']
        dmat = refs['d_matrix']
        predictor.for_reflection_table(refs, UB, s0, dmat)
      # scan static
      else:
        predictor = sc(e)
        UB = e.crystal.get_A()
        predictor.for_reflection_table(refs, UB)

      # write predictions back to overall reflections
      reflections.set_selected(sel, refs)

    return reflections

