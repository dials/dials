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

from __future__ import division

from math import pi
from scitbx import matrix
from dials.algorithms.spot_prediction import RayPredictor

from dials.model.data import Reflection, ReflectionList

from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
from dials.algorithms.spot_prediction import StillsReflectionPredictor


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
    self.update()

  def update(self):
    """Build RayPredictor objects for the current geometry of each Experiment"""

    self._ray_predictors = [RayPredictor(
                              e.beam.get_s0(),
                              e.goniometer.get_rotation_axis(),
                              self._sweep_range) for e in self._experiments]
    self._UBs = [e.crystal.get_U() * e.crystal.get_B() for e in self._experiments]

  def predict(self, hkl, experiment_id=0, UB=None):
    """
    Solve the prediction formula for the reflecting angle phi.

    If UB is given, override the contained crystal model. This is
    for use in refinement with time-varying crystal parameters
    """

    UB_ = UB if UB else self._UBs[experiment_id]

    return self._ray_predictors[experiment_id](hkl, UB_)

class StillsRayPredictor(object):
  """
  Predict for a relp based on the current states of models of the
  experimental geometry. Here we assume the crystal UB already puts hkl
  in reflecting position, so no rotation is required.

  Generally, this assumption is not true: most relps are not touching the
  Ewald sphere on a still image, but are observed due to their finite
  mosaicity and the finite width of the Ewald sphere. Rather than employing
  a complicated polychromatic and mosaic model, here we take a naive
  approach, which is likely to introduce (small?) errors in the direction of
  predicted rays.

  """

  def __init__(self, experiments):
    """Construct by linking to instances of experimental model classes"""

    self._experiments = experiments
    self.update()

  def update(self):
    """Cache s0 and UB"""

    self._s0s = [matrix.col(e.beam.get_s0()) for e in self._experiments]
    self._s0_lengths = [s0.length() for s0 in self._s0s]
    self._UBs = [e.crystal.get_U() * e.crystal.get_B() for e in self._experiments]

  def predict(self, hkl, experiment_id=0):
    """Predict for hkl under the assumption it is in reflecting position"""

    s0 = self._s0s[experiment_id]
    s0_length = self._s0_lengths[experiment_id]
    r = self._UBs[experiment_id] * matrix.col(hkl)
    s1 = (s0 + r).normalize() * s0_length

    # create the Reflections and set properties. The relp is
    # neither entering nor exiting the Ewald sphere, but we need
    # both to ensure we match whatever this flag is set to for
    # the observation.
    ref1, ref2 = Reflection(hkl), Reflection(hkl)
    ref1.beam_vector, ref2.beam_vector = s1, s1
    ref1.rotation_angle, ref2.rotation_angle = 0.0, 0.0
    ref1.entering, ref2.entering = True, False

    rl = ReflectionList(2)
    rl[0] = ref1
    rl[1] = ref2
    return rl

class ExperimentsPredictor(object):
  """
  Predict for relps based on the current states of models of the experimental
  geometry. This version manages multiple experiments, selecting the correct
  predictor in each case.
  """

  def __init__(self, experiments):
    """Construct by linking to instances of experimental model classes"""

    self._experiments = experiments
    self.update()

  def update(self):
    """Build predictor objects for the current geometry of each Experiment"""

    sc = ScanStaticReflectionPredictor
    st = StillsReflectionPredictor
    self._predictors = [sc(e) if e.goniometer else st(e) \
                        for e in self._experiments]
    self._UBs = [e.crystal.get_U() * e.crystal.get_B() for e in self._experiments]

  def predict(self, reflections):
    """
    Predict for all reflections
    """

    for iexp, e in enumerate(self._experiments):

      # select the reflections for this experiment only
      sel = reflections['id'] == iexp
      refs = reflections.select(sel)

      # determine whether to try scan-varying prediction
      if refs.has_key('ub_matrix'):
        UBs = refs['ub_matrix']
        # predict and assign in place
        self._predictors[iexp].for_reflection_table(refs, UBs)
      else:
        # predict and assign in place
        self._predictors[iexp].for_reflection_table(refs, self._UBs[iexp])

      # write predictions back to overall reflections
      reflections.set_selected(sel, refs)

    return reflections

