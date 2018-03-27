#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

""" Test the situation that led to https://github.com/dials/dials/issues/423.
In that case instantiating a Refiner for an experiment list with an I23
detector model caused the panel origins to move before any refinement took
place. This occured because for the input experiments.json the root frame for
the hierarchical detector is on source side of the laboratory frame origin, not
on the detector side. Prior to the fix this resulted in incorrect calculation
of the offsets of all panels from the root frame.
"""

from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory

class Test(object):

  def __init__(self):

    dials_regression = libtbx.env.find_in_repositories(
      relative_path="dials_regression",
      test=os.path.isdir)

    data_dir = os.path.join(dials_regression, "refinement_test_data",
        "dials-423")
    exp_file = os.path.join(data_dir, 'experiments.json')
    ref_file = os.path.join(data_dir, 'subset.pickle')

    self._reflections = flex.reflection_table.from_pickle(ref_file)
    self._experiments = ExperimentListFactory.from_json_file(exp_file,
        check_format=False)

  def run(self):
    """Test that the detector remains similar after refiner construction"""

    from dials.algorithms.refinement.refiner import phil_scope
    params = phil_scope.fetch(source=phil.parse('')).extract()

    # disable outlier rejection for speed of refiner construction
    params.refinement.reflections.outlier.algorithm='null'

    refiner = RefinerFactory.from_parameters_data_experiments(params,
        self._reflections, self._experiments)

    d1 = self._experiments[0].detector
    d2 = refiner.get_experiments()[0].detector

    assert d1.is_similar_to(d2)
    return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  tst = Test()
  tst.run()

if __name__ == '__main__':
  run()
