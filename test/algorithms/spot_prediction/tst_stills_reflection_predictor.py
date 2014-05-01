#!/usr/bin/env python
#
#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from scitbx import matrix

class Test(object):

  def __init__(self):

    from dials.model.experiment import beam_factory
    from dials.model.experiment import goniometer_factory
    from dials.model.experiment import detector_factory

    from cctbx.crystal.crystal_model import crystal_model

    self.beam = beam_factory.make_beam(unit_s0 = matrix.col((0, 0, 1)),
                                       wavelength = 1.0)
    #self.goniometer = goniometer_factory.known_axis(matrix.col((1, 0, 0)))

    dir1 = matrix.col((1, 0, 0))
    dir2 = matrix.col((0, -1, 0))
    n = matrix.col((0, 0, 1))
    centre = matrix.col((0, 0, 200))
    npx_fast = npx_slow = 1000
    pix_size = 0.2
    origin = centre - (0.5 * npx_fast * pix_size * dir1 +
                       0.5 * npx_slow * pix_size * dir2)
    self.detector = detector_factory.make_detector("PAD",
                        dir1, dir2, origin,
                        (pix_size, pix_size),
                        (npx_fast, npx_slow),
                        (0, 1.e6))

    a = matrix.col((1, 0, 0))
    b = matrix.col((0, 1, 0))
    c = matrix.col((0, 0, 1))
    self.crystal = crystal_model(a, b, c, space_group_symbol = "P 1")

    from dials.model.experiment.experiment_list import ExperimentList, Experiment
    self.experiments = ExperimentList()
    self.experiments.append(Experiment(
      beam=self.beam, detector=self.detector, goniometer=None,
      scan=None, crystal=self.crystal, imageset=None))

    from dials.array_family import flex
    table = flex.reflection_table()
    table['flags'] = flex.size_t([0])
    table['id']    = flex.size_t([0])
    table['panel'] = flex.size_t([0])

    # Predicted properties
    table['miller_index'] = flex.miller_index([(1,0,0)])
    table['entering']     = flex.bool([True])
    table['s1']           = flex.vec3_double(1)
    table['xyzcal.mm']    = flex.vec3_double(1)
    table['xyzcal.px']    = flex.vec3_double(1)

    self.reflections = table
    print self.reflections[0]
    for k in self.reflections.keys(): print k


  def run(self):

    from dials.algorithms.spot_prediction import StillsReflectionPredictor

    predictor = StillsReflectionPredictor(self.experiments[0])

    tmp = predictor.for_reflection_table(self.reflections)


if __name__ == '__main__':

  test = Test()
  test.run()
