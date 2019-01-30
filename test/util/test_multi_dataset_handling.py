"""
Tests for dials.util.multi_dataset_handling functions
"""

import pytest
from dials.util import Sorry
from dxtbx.model import Experiment, ExperimentList
from dials.array_family import flex
from dials.util.multi_dataset_handling import assign_unique_identifiers,\
  parse_multiple_datasets, select_datasets_on_ids

def empty_explist_3():
  """Make a list of three empty experiments"""
  experiments = ExperimentList()
  experiments.append(Experiment())
  experiments.append(Experiment())
  experiments.append(Experiment())
  return experiments

def reflection_list_3():
  """Make a list of three reflection tables"""
  rt1 = flex.reflection_table()
  rt1['id'] = flex.int([0, 0, 0])
  rt2 = flex.reflection_table()
  rt2['id'] = flex.int([1, 1])
  rt3 = flex.reflection_table()
  rt3['id'] = flex.int([4, 4])
  reflections = [rt1, rt2, rt3]
  return reflections

def multi_table_3():
  r = flex.reflection_table()
  r['id'] = flex.int([0, 1, 4])
  return [r]

def test_select_datasets_on_ids():
  """Test the select_datasets_on_ids function."""
  experiments = empty_explist_3()
  reflections = reflection_list_3()
  reflections[0].experiment_identifiers()[0] = '0'
  experiments[0].identifier = '0'
  reflections[1].experiment_identifiers()[1] = '2'
  experiments[1].identifier = '2'
  reflections[2].experiment_identifiers()[4] = '4'
  experiments[2].identifier = '4'
  use_datasets = ['0', '2']
  experiments, refl = select_datasets_on_ids(experiments, reflections,
    use_datasets=use_datasets)
  assert len(experiments) == 2
  assert len(refl) == 2
  assert list(experiments.identifiers()) == ['0', '2']

  experiments = empty_explist_3()
  experiments[0].identifier = '0'
  experiments[1].identifier = '2'
  experiments[2].identifier = '4'
  exclude_datasets = ['0']
  experiments, refl = select_datasets_on_ids(experiments, reflections,
    exclude_datasets=exclude_datasets)
  assert len(refl) == 2
  assert list(experiments.identifiers()) == ['2', '4']
  assert len(experiments) == 2

  with pytest.raises(Sorry):
    exclude_datasets = ['0']
    use_datasets = ['2', '4']
    experiments, refl = select_datasets_on_ids(experiments,
      reflections, use_datasets=use_datasets,
      exclude_datasets=exclude_datasets)

  with pytest.raises(Sorry):
    exclude_datasets = ['1']
    experiments, refl = select_datasets_on_ids(experiments,
      reflections, exclude_datasets=exclude_datasets)

  # test for return if no identifiers specified
  exp, refl = select_datasets_on_ids(experiments, reflections)
  assert exp is experiments
  assert refl is reflections

  # test for Sorry raise if not all identifiers set
  experiments = empty_explist_3()
  experiments[0].identifier = '0'
  experiments[1].identifier = '2'
  with pytest.raises(Sorry):
    exp, refl = select_datasets_on_ids(experiments, reflections, use_datasets=['2'])
  # test for Sorry raise if bad input
  experiments[2].identifier = '4'
  with pytest.raises(Sorry):
    exp, refl = select_datasets_on_ids(experiments, reflections, use_datasets=['3'])

  # test correct handling with multi-dataset table
  reflections = flex.reflection_table()
  reflections['id'] = flex.int([0, 1, 2])
  reflections.experiment_identifiers()[0] = '0'
  reflections.experiment_identifiers()[1] = '2'
  reflections.experiment_identifiers()[2] = '4'
  exp, refl = select_datasets_on_ids(experiments, [reflections], exclude_datasets=['2'])
  assert list(refl[0].experiment_identifiers().values()) == ['0', '4']
  assert list(refl[0]['id']) == [0, 2]

def test_assign_unique_identifiers():
  """Test the assignment of unique identifiers"""
  # Test case where none are set but refl table ids are - use refl ids
  experiments = empty_explist_3()
  reflections = reflection_list_3()
  assert list(experiments.identifiers()) == ['', '', '']
  exp, rts = assign_unique_identifiers(experiments, reflections)
  expected_identifiers = ['0', '1', '2']
  expected_ids = [0, 1, 2]
  # Check that identifiers are set in experiments and reflection table.
  assert (list(exp.identifiers())) == expected_identifiers
  for i, refl in enumerate(rts):
    assert refl.experiment_identifiers()[i] == expected_identifiers[i]
    assert list(set(refl['id'])) == [i]
    assert i == expected_ids[i]

  # Test case where none are set but refl table ids have duplicates
  experiments = empty_explist_3()
  reflections = reflection_list_3()
  reflections[2]['id'] = flex.int([0, 0])
  exp, rts = assign_unique_identifiers(experiments, reflections)
  expected_identifiers = ['0', '1', '2']
  # Check that identifiers are set in experiments and reflection table.
  assert (list(exp.identifiers())) == expected_identifiers
  expected_ids = [0, 1, 2]
  for i, refl in enumerate(rts):
    assert refl.experiment_identifiers()[i] == expected_identifiers[i]
    assert list(set(refl['id'])) == [i]
    assert i == expected_ids[i]

  # Test case where identifiers are already set.
  experiments = empty_explist_3()
  experiments[0].identifier = '0'
  experiments[1].identifier = '4'
  experiments[2].identifier = '2'
  reflections = reflection_list_3()
  reflections[1].experiment_identifiers()[0] = '5'
  # should raise an assertion error for inconsistent identifiers
  with pytest.raises(ValueError):
    exp, rts = assign_unique_identifiers(experiments, reflections)

  #test cases where all set, whether reflection table is split or not
  experiments = empty_explist_3()
  experiments[0].identifier = '0'
  experiments[1].identifier = '4'
  experiments[2].identifier = '2'
  reflections = reflection_list_3()
  reflections[0].experiment_identifiers()[0] = '0'
  reflections[1].experiment_identifiers()[1] = '4'
  reflections[2].experiment_identifiers()[4] = '2'

  experiments_multi = empty_explist_3()
  reflections_multi = multi_table_3()
  experiments_multi[0].identifier = '0'
  experiments_multi[1].identifier = '4'
  experiments_multi[2].identifier = '2'
  reflections_multi[0].experiment_identifiers()[0] = '0'
  reflections_multi[0].experiment_identifiers()[1] = '4'
  reflections_multi[0].experiment_identifiers()[4] = '2'
  #should pass experiments back if identifiers all already set
  for exper, refl in zip([experiments, experiments_multi], [reflections, reflections_multi]):
    exp, rts = assign_unique_identifiers(exper, refl)
    expected_identifiers = ['0', '4', '2']
    # Check that identifiers are set in experiments and reflection table.
    assert exp is exper
    assert list(exp.identifiers()) == expected_identifiers
    expected_ids = [0, 1, 4]
    for i, refl in enumerate(rts):
      id_ = refl['id'][0]
      assert refl.experiment_identifiers()[id_] == expected_identifiers[i]
      assert list(set(refl['id'])) == [id_]
      assert id_ == expected_ids[i]

  # Test for correct Sorry raise if unqeual experiments and reflections
  del reflections_multi[0].experiment_identifiers()[4]
  del reflections_multi[0]['id'][2]
  with pytest.raises(Sorry):
    exp, rts = assign_unique_identifiers(experiments_multi, reflections_multi)

  # Now test that if some are set, these are maintained and unique ids are
  # set for the rest
  experiments = empty_explist_3()
  experiments[0].identifier = '1'
  reflections = reflection_list_3()
  reflections[0].experiment_identifiers()[0] = '1'
  exp, rts = assign_unique_identifiers(experiments, reflections)
  expected_identifiers = ['1', '0', '2']
  assert list(exp.identifiers()) == expected_identifiers
  expected_ids = [0, 1, 2]
  for i, refl in enumerate(rts):
    id_ = refl['id'][0]
    assert refl.experiment_identifiers()[id_] == expected_identifiers[i]
    assert list(set(refl['id'])) == [id_]
    assert id_ == expected_ids[i]

  # Test the case where identifiers are specified
  experiments = empty_explist_3()
  reflections = reflection_list_3()
  reflections[0].experiment_identifiers()[0] = '5'
  reflections[1].experiment_identifiers()[1] = '6'
  reflections[1].experiment_identifiers()[4] = '7'
  input_identifiers = ['0', '1', '10']
  exp, rts = assign_unique_identifiers(experiments, reflections, input_identifiers)
  assert list(exp.identifiers()) == input_identifiers
  assert rts[0].experiment_identifiers()[0] == '0'
  assert rts[0]['id'][0] == 0
  assert rts[1].experiment_identifiers()[1] == '1'
  assert rts[1]['id'][0] == 1
  assert rts[2].experiment_identifiers()[2] == '10'
  assert rts[2]['id'][0] == 2

  # Test raises Sorry when wrong number of identifiers given
  with pytest.raises(Sorry):
    exp, rts = assign_unique_identifiers(experiments, reflections, ['0', '1'])

def test_parse_multiple_datasets():
  """Test the namesake function. This expects a list of reflection tables, and
  selects on the column named 'id'."""
  rt1 = flex.reflection_table()
  rt1['id'] = flex.int([0, 0, 0])
  rt1.experiment_identifiers()[0] = '0'
  rt2 = flex.reflection_table()
  rt2['id'] = flex.int([2, 2, 4, 4])
  rt2.experiment_identifiers()[2] = '2'
  rt2.experiment_identifiers()[4] = '4'
  single_tables = parse_multiple_datasets([rt2])
  assert len(single_tables) == 2
  assert list(set(single_tables[0]['id'])) == [2]
  assert list(set(single_tables[1]['id'])) == [4]
  single_tables = parse_multiple_datasets([rt1, rt2])
  assert list(set(single_tables[0]['id'])) == [0]
  assert list(set(single_tables[1]['id'])) == [2]
  assert list(set(single_tables[2]['id'])) == [4]
  assert len(single_tables) == 3
  single_tables = parse_multiple_datasets([rt1])
  assert len(single_tables) == 1
  assert list(set(single_tables[0]['id'])) == [0]

  # if a duplicate id is given, then this should be detected and new ids
  # determined for all datasets.
  rt3 = flex.reflection_table()
  rt3['id'] = flex.int([2, 2])
  rt3.experiment_identifiers()[2] = '5'
  single_tables = parse_multiple_datasets([rt1, rt2, rt3])
  assert len(single_tables) == 4
  assert list(set(single_tables[0]['id'])) == [0]
  assert single_tables[0].experiment_identifiers()[0] == '0'
  assert list(set(single_tables[1]['id'])) == [1]
  assert single_tables[1].experiment_identifiers()[1] == '2'
  assert list(set(single_tables[2]['id'])) == [2]
  assert single_tables[2].experiment_identifiers()[2] == '4'
  assert list(set(single_tables[3]['id'])) == [3]
  assert single_tables[3].experiment_identifiers()[3] == '5'
