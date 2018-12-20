from copy import deepcopy
import pytest
from mock import Mock
from dials.algorithms.scaling.simple_Ih_table import simple_Ih_table, Ih_table_block
from dials.array_family import flex
from cctbx.sgtbx import space_group
from scitbx import sparse

@pytest.fixture()
def large_reflection_table():
  """Create a reflection table."""
  return generate_refl_1()

@pytest.fixture()
def small_reflection_table():
  """Create a small reflection table."""
  return generate_refl_2()

@pytest.fixture(scope='module')
def test_sg():
  """Create a space group object."""
  return space_group("C 2y")

def generate_refl_1():
  """Generate a test reflection table. Note tha the variance values are chosen
  as the 'True' Ih_values, which would be found if unity weights were chosen
  in this example."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([100.0, 100.0, 80.0, 60.0, 30.0,
    40.0, 60.0])
  reflections['inverse_scale_factor'] = flex.double(7, 1.0)
  reflections['variance'] = flex.double([90.0, 100.0, 90.0, 60.0, 30.0,
    50.0, 50.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (-1, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, -2), (0, 0, 2)])
  reflections.set_flags(flex.bool(7, True), reflections.flags.integrated)
  return reflections

def generate_refl_2():
  """Generate another test reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([60.0, 30.0, 10.0, 30.0])
  reflections['variance'] = flex.double([60.0, 30.0, 10.0, 30.0])
  reflections['inverse_scale_factor'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 4, 0),
    (10, 0, 0), (0, 4, 0)])
  reflections.set_flags(flex.bool(4, True), reflections.flags.integrated)
  return reflections

def test_IhTable_split_into_blocks(large_reflection_table,
  small_reflection_table, test_sg):
  """Test that the Ih_table datastructure correctly organises the data
  from two reflection tables into two IhTableBlocks."""

  sel1 = flex.bool(7, True)
  sel1[6] = False
  sel2 = flex.bool(4, True)
  sel2[1] = False

  Ih_table = simple_Ih_table([(large_reflection_table, sel1), (small_reflection_table, sel2)],
    test_sg, nblocks=2)

  # reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
  #  (-1, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, -2), (0, 0, 2)])
  # reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 4, 0),
  #  (10, 0, 0), (0, 4, 0)])

  assert Ih_table.n_datasets == 2
  assert Ih_table.nblocks == 2
  block_list = Ih_table.Ih_table_blocks
  assert list(block_list[0].Ih_table['asu_miller_index']) == [(0, 0, 1),
    (0, 0, 2), (0, 2, 0)]
  assert list(block_list[1].Ih_table['asu_miller_index']) == [(0, 4, 0),
    (1, 0, 0), (1, 0, 0), (0, 4, 0), (1, 0, 0), (10, 0, 0)]
  assert list(block_list[0].block_selections[0]) == [1, 5, 3]
  assert list(block_list[0].block_selections[1]) == []
  assert list(block_list[1].block_selections[0]) == [4, 0, 2]
  assert list(block_list[1].block_selections[1]) == [3, 0, 2]

  expected_h_idx_matrix = sparse.matrix(4, 3)
  expected_h_idx_matrix[0, 0] = 1
  expected_h_idx_matrix[1, 1] = 1
  expected_h_idx_matrix[2, 2] = 1
  assert block_list[0].h_index_matrix == expected_h_idx_matrix
  expected_h_idx_matrix = sparse.matrix(7, 3)
  expected_h_idx_matrix[0, 0] = 1
  expected_h_idx_matrix[1, 1] = 1
  expected_h_idx_matrix[2, 1] = 1
  expected_h_idx_matrix[3, 0] = 1
  expected_h_idx_matrix[4, 1] = 1
  expected_h_idx_matrix[5, 2] = 1
  assert block_list[1].h_index_matrix == expected_h_idx_matrix
