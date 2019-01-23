# coding: utf-8
"""
Testing for other refiner/refinement small utilities that can be tested
alone.
"""

from __future__ import absolute_import, division, print_function

from mock import Mock, patch
from dials.algorithms.refinement.refiner import _copy_experiments_for_refining

def test_that_the_experiment_reduced_copier_works_as_intended():
  sample = Mock()
  with patch('dials.algorithms.refinement.refiner.ExperimentList', new=list):
    dupe = _copy_experiments_for_refining(sample)[0]
  assert dupe is not sample
  # Make sure that anything refined is unique
  for att in ["beam", "goniometer", "detector", "crystal"]:
    assert getattr(sample, att) is not getattr(dupe, att)
  # Anything read-only should be untouched
  for att in ["scan", "profile", "imageset", "scaling_model"]:
    assert getattr(sample, att) is getattr(dupe, att)
