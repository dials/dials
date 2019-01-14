"""
Tests for the functions in dials.util.options
"""
from mock import Mock
from dials.util.options import flatten_reflections
from dials.array_family import flex

def mock_reflection_file_object(id_=0, identifier=True):
  """Create a mock reflection_file_object."""
  fileobj = Mock()
  r = flex.reflection_table()
  r['id'] = flex.int([-1, id_, id_])
  if identifier:
    r.experiment_identifiers()[id_] = str(id_)
  fileobj.data = r
  return fileobj

def mock_two_reflection_file_object():
  """Create a mock reflection_file_object with two datasets."""
  fileobj = Mock()
  r = flex.reflection_table()
  r['id'] = flex.int([-1, 0, 0, 2, 2])
  r.experiment_identifiers()[0] = str(0)
  r.experiment_identifiers()[2] = str(2)
  fileobj.data = r
  return fileobj

def test_flatten_experiments_updating_id_values():
  """Test the correct handling of duplicate table id values.

  Note that this function does not have the ability to update the
  experiment string identifier, only ensure that the table id values
  do not clash (it is not possible even to load multiple experiments
  with the same identifier).
  """
  # Test the case of two single reflection tables.
  file_list = [mock_reflection_file_object(id_=0),
    mock_reflection_file_object(id_=0)]
  rs = flatten_reflections(file_list)
  assert rs[0] is file_list[0].data
  assert list(rs[0]['id']) == [-1, 0, 0]
  assert list(rs[0].experiment_identifiers().keys()) == [0]
  assert list(rs[0].experiment_identifiers().values()) == ['0']
  assert rs[1] is file_list[1].data
  assert list(rs[1]['id']) == [-1, 1, 1]
  assert list(rs[1].experiment_identifiers().keys()) == [1]
  assert list(rs[1].experiment_identifiers().values()) == ['0']

  # Now test the case where one reflection table contains two experiments
  file_list = [mock_two_reflection_file_object(),
    mock_reflection_file_object(id_=0)]
  rs = flatten_reflections(file_list)
  assert rs[0] is file_list[0].data
  assert list(rs[0]['id']) == [-1, 0, 0, 1, 1]
  assert list(rs[0].experiment_identifiers().keys()) == [0, 1]
  assert list(rs[0].experiment_identifiers().values()) == ['0', '2']
  assert rs[1] is file_list[1].data
  assert list(rs[1]['id']) == [-1, 2, 2]
  assert list(rs[1].experiment_identifiers().keys()) == [2]
  assert list(rs[1].experiment_identifiers().values()) == ['0']
