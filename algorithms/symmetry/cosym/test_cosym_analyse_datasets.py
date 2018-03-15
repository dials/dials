from __future__ import absolute_import, division, print_function

import pytest

from cctbx import sgtbx

from dials.algorithms.symmetry.cosym.generate_test_data import generate_test_data
from dials.algorithms.symmetry.cosym import phil_scope
from dials.algorithms.symmetry.cosym import analyse_datasets


@pytest.mark.parametrize('space_group', ['P2', 'P3', 'I23'])
def test_cosym_analyse_datasets(space_group):
  datasets, expected_reindexing_ops = generate_test_data(
    space_group=sgtbx.space_group_info(symbol=space_group).group())
  expected_space_group = sgtbx.space_group_info(symbol=space_group).group()

  params = phil_scope.extract()
  params.cluster.agglomerative.n_clusters = len(expected_reindexing_ops)

  result = analyse_datasets(datasets, params)

  space_groups = {}
  reindexing_ops = {}
  for dataset_id in result.reindexing_ops.iterkeys():
    if 0 in result.reindexing_ops[dataset_id]:
      cb_op = result.reindexing_ops[dataset_id][0]
      reindexing_ops.setdefault(cb_op, set())
      reindexing_ops[cb_op].add(dataset_id)
    if dataset_id in result.space_groups:
      space_groups.setdefault(result.space_groups[dataset_id], set())
      space_groups[result.space_groups[dataset_id]].add(dataset_id)

  assert len(reindexing_ops) == len(expected_reindexing_ops)
  assert sorted(reindexing_ops.keys()) == sorted(expected_reindexing_ops.keys())

  for ridx_set in reindexing_ops.values():
    for expected_set in expected_reindexing_ops.values():
      assert (
        #(len(ridx_set.symmetric_difference(expected_set)) <= 2) or
        #(len(ridx_set.intersection(expected_set)) <= 2))
        (len(ridx_set.symmetric_difference(expected_set)) == 0) or
        (len(ridx_set.intersection(expected_set)) == 0))


  #assert expected_reindexing_ops == reindexing_ops

  #for op in expected_reindexing_ops:
    #assert op in reindexing_ops, op
    #assert expected_reindexing_ops
