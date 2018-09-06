'''
This code tests for Ih_table and joint_Ih_table data structures.
This also provides a test for the scaler, which must be successfully
initialised in order to provide input for the Ih_table.
'''
from copy import deepcopy
import pytest
from mock import Mock
from dials.algorithms.scaling.Ih_table import IhTable, IhTableBlock
from dials.array_family import flex
from cctbx.sgtbx import space_group
from scitbx import sparse

@pytest.fixture()
def reflection_table_for_block():
  """Create a reflection table."""
  return generated_refl_for_splitting_1()

@pytest.fixture()
def large_reflection_table():
  """Create a reflection table."""
  return generate_refl_1()

def generated_refl_for_splitting_1():
  """"Make a reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
  reflections['variance'] = flex.double(6, 1.0)
  reflections['inverse_scale_factor'] = flex.double(6, 1.0)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (2, 0, 0), (0, 0, 1),
    (2, 2, 2), (1, 0, 0), (2, 0, 0)])
  reflections['d'] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6, 2.5])
  reflections['partiality'] = flex.double(6, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 8.0), (0.0, 0.0, 10.0), (0.0, 0.0, 12.0),
    (0.0, 0.0, 15.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0), (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool(6, True), reflections.flags.integrated)
  reflections.set_flags(flex.bool(6, False), reflections.flags.bad_for_scaling)
  return [reflections]

@pytest.fixture()
def refl_for_target_Ih_table_1():
  """A target reflection table for 'generated_refl_for_splitting_1'"""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([10.0, 20.0])
  reflections['variance'] = flex.double(2, 1.0)
  reflections['inverse_scale_factor'] = flex.double(2, 1.0)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (2, 0, 0)])
  reflections.set_flags(flex.bool(2, True), reflections.flags.integrated)
  reflections.set_flags(flex.bool(2, False), reflections.flags.bad_for_scaling)
  return reflections

@pytest.fixture()
def refl_for_target_Ih_table_2():
  """A target reflection table for 'generated_refl_for_splitting_1'"""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([50.0, 10.0])
  reflections['variance'] = flex.double(2, 1.0)
  reflections['inverse_scale_factor'] = flex.double(2, 1.0)
  reflections['miller_index'] = flex.miller_index([(2, 2, 2), (2, 0, 0)])
  reflections.set_flags(flex.bool(2, True), reflections.flags.integrated)
  reflections.set_flags(flex.bool(2, False), reflections.flags.bad_for_scaling)
  return reflections

def generated_refl_for_splitting_2():
  """Make a reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([7.0, 8.0, 9.0, 10.0, 11.0])
  reflections['variance'] = flex.double(5, 1.0)
  reflections['inverse_scale_factor'] = flex.double(5, 1.0)
  reflections['miller_index'] = flex.miller_index([(2, 2, 2), (2, 0, 0), (0, 0, 1),
    (2, 2, 2), (1, 0, 0)])
  reflections['d'] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6])
  reflections['partiality'] = flex.double(5, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 8.0), (0.0, 0.0, 10.0), (0.0, 0.0, 12.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool(5, True), reflections.flags.integrated)
  reflections.set_flags(flex.bool(5, False), reflections.flags.bad_for_scaling)
  return [reflections]

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

def mock_error_model():
  """Mock error model."""
  em = Mock()
  em.update_variances.return_value = flex.double([1.0, 2.0, 3.0])
  return em

def mock_error_model_2():
  """Mock error model."""
  em = Mock()
  em.update_variances.return_value = flex.double([1.0, 2.0, 3.0, 1.0, 2.0, 3.0])
  return em

def test_IhTableBlock(reflection_table_for_block, test_sg):
  """A test for the IhTableBlock class, the main datastructure
  for managing a reflection table with an indexing matrix to allow summation
  over groups of equivalent reflections."""
  # The block should maintain the initial order of the reflection table
  rt = reflection_table_for_block[0]
  Ihtableblock = IhTableBlock(rt, test_sg)

  assert Ihtableblock.space_group is test_sg

  # Test that the asu miller index was set and the correct h_idx matrix was set.
  assert list(Ihtableblock.asu_miller_index) == [(1, 0, 0), (2, 0, 0),
    (0, 0, 1), (2, 2, 2), (1, 0, 0), (2, 0, 0)]
  h_idx = Ihtableblock.h_index_matrix
  assert h_idx[2, 0] == 1
  assert h_idx[0, 1] == 1
  assert h_idx[4, 1] == 1
  assert h_idx[1, 2] == 1
  assert h_idx[5, 2] == 1
  assert h_idx[3, 3] == 1
  assert h_idx.n_cols == 4
  assert h_idx.n_rows == 6
  assert h_idx.non_zeroes == 6

  h_expand = Ihtableblock.h_expand_matrix
  assert h_expand[0, 2] == 1
  assert h_expand[1, 0] == 1
  assert h_expand[1, 4] == 1
  assert h_expand[2, 1] == 1
  assert h_expand[2, 5] == 1
  assert h_expand[3, 3] == 1
  assert h_expand.n_cols == 6
  assert h_expand.n_rows == 4
  assert h_idx.non_zeroes == 6

  assert Ihtableblock.derivatives is None
  assert Ihtableblock.nonzero_weights is None
  assert Ihtableblock.n_h is None

  # Test other expected properties
  assert Ihtableblock.size == 6
  assert list(Ihtableblock.weights) == [1.0] * 6
  assert list(Ihtableblock.variances) == list(rt['variance'])
  assert list(Ihtableblock.intensities) == list(rt['intensity'])
  assert list(Ihtableblock.inverse_scale_factors) == [1.0] * 6
  assert list(Ihtableblock.Ih_values) == [3.0, 4.0, 3.0, 4.0, 3.0, 4.0]
  assert Ihtableblock.Ih_table
  assert list(Ihtableblock.miller_index) == list(rt['miller_index'])

  # Test calc n_h
  Ihtableblock.calc_nh()
  assert list(Ihtableblock.n_h) == [2, 2, 1, 1, 2, 2]

  # Test calc Ih - set new scale factors and check Ih values are
  # correctly calculated.
  Ihtableblock.inverse_scale_factors = flex.double([6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
  Ihtableblock.calc_Ih()
  assert list(Ihtableblock.Ih_values) == [16.0/40.0, 16.0/26.0, 3.0/4.0,
    4.0/3.0, 16.0/40.0, 16.0/26.0]

  # Test set_Ih_values_to_target, do it twice with different targets just in case
  target_refl_1 = refl_for_target_Ih_table_1()
  target_Ih_table_1 = IhTableBlock(target_refl_1, test_sg)
  target_refl_2 = refl_for_target_Ih_table_2()
  target_Ih_table_2 = IhTableBlock(target_refl_2, test_sg)

  Ihtableblock.set_Ih_values_to_target(target_Ih_table_1)
  assert list(Ihtableblock.Ih_values) == [10.0, 20.0, 0.0, 0.0, 10.0, 20.0]
  Ihtableblock.set_Ih_values_to_target(target_Ih_table_2)
  assert list(Ihtableblock.Ih_values) == [0.0, 10.0, 0.0, 50.0, 0.0, 10.0]

  # Now test updating certain properties and that appropriate errors are raised
  with pytest.raises(AssertionError):
    Ihtableblock.inverse_scale_factors = flex.double([1.0, 2.0])
  Ihtableblock.inverse_scale_factors = flex.double(6, 2.0)
  assert list(Ihtableblock.inverse_scale_factors) == [2.0] * 6

  with pytest.raises(AssertionError):
    Ihtableblock.weights = flex.double([1.0, 2.0])
  Ihtableblock.weights = flex.double(6, 3.0)
  assert list(Ihtableblock.weights) == [3.0] * 6

  with pytest.raises(AssertionError):
    Ihtableblock.h_index_matrix = sparse.matrix(1, 1)
  with pytest.raises(AssertionError):
    Ihtableblock.h_expand_matrix = sparse.matrix(1, 1)

  with pytest.raises(AssertionError):
    Ihtableblock.nonzero_weights = flex.int([1])
  Ihtableblock.nonzero_weights = flex.int([5, 6, 7, 8, 9, 10])
  assert list(Ihtableblock.nonzero_weights) == [5, 6, 7, 8, 9, 10]

  # Test select
  Ihtableblock.select(flex.bool([True, False, False, True, True, False]))
  assert list(Ihtableblock.asu_miller_index) == [(1, 0, 0), (2, 2, 2), (1, 0, 0)]
  h_idx = Ihtableblock.h_index_matrix
  assert h_idx[0, 0] == 1
  assert h_idx[1, 1] == 1
  assert h_idx[2, 0] == 1
  assert h_idx.n_cols == 2
  assert h_idx.n_rows == 3
  assert h_idx.non_zeroes == 3

  h_expand = Ihtableblock.h_expand_matrix
  assert h_expand[0, 0] == 1
  assert h_expand[1, 1] == 1
  assert h_expand[0, 2] == 1
  assert h_expand.n_cols == 3
  assert h_expand.n_rows == 2
  assert h_expand.non_zeroes == 3

  assert list(Ihtableblock.nonzero_weights) == [5, 8, 9]

  #Test update error model
  Ihtableblock.update_error_model(mock_error_model())
  assert list(Ihtableblock.weights) == [1.0, 1.0/2.0, 1.0/3.0]

def test_apply_iterative_weighting(reflection_table_for_block, test_sg):
  """Test the setting of iterative weights."""

  Ihtableblock = IhTableBlock(reflection_table_for_block[0], test_sg,
    weighting_scheme='GM')

  # After construction, weights should initially be set to inverse variances,
  # to allow the first calculation of Ih (by least-squares approach).
  assert list(Ihtableblock.weights) == [1.0] * 6

  # Now update weights
  Ihtableblock.update_weights()
  assert list(Ihtableblock.weights) != [1.0] * 6 # Check changed, or that
  # new weights are not by chance identical to old weights.
  gIh = Ihtableblock.inverse_scale_factors * Ihtableblock.Ih_values
  t = (Ihtableblock.intensities - gIh) / gIh
  assert list(Ihtableblock.weights) == list(1.0/(1.0 + t**2)**2)


  Ih_table = IhTable([(reflection_table_for_block[0], None)], test_sg,
    weighting_scheme='GM')
  block = Ih_table.blocked_data_list[0]
  assert list(block.weights) == [1.0] * 6
  Ih_table.update_weights()
  gIh = block.inverse_scale_factors * block.Ih_values
  t = (block.intensities - gIh) / gIh
  assert list(block.weights) == list(1.0/(1.0 + t**2)**2)

def test_IhTable(reflection_table_for_block, test_sg):
  """Test for Ih_table datastructure. Upon initialisation, Ih_table should set
  unity scale factors and calculate Ih_values. It should also create the
  a h_index_matrix.
  This will be similar to the first test, but this time the data is sorted
  by asu miller index."""
  large_reflection_table = reflection_table_for_block[0]

  # Test initialisation fails without correct columns
  for key in ['intensity', 'variance', 'inverse_scale_factor', 'miller_index']:
    col_copy = deepcopy(large_reflection_table[key])
    del large_reflection_table[key]
    with pytest.raises(AssertionError):
      Ih_table = IhTable([(large_reflection_table, None)], test_sg,
        n_blocks=1, weighting_scheme='unity')
    large_reflection_table[key] = col_copy

  Ih_table = IhTable([(large_reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')

  assert Ih_table.id_ == "IhTable"
  assert Ih_table.n_datasets == 1

  # Tests calc_Ih, assign_h_matrices, interface
  block = Ih_table.blocked_data_list[0]
  assert block.size == 6
  assert list(block.inverse_scale_factors) == [1.0] * 6
  assert list(block.weights) == [1.0] * 6
  assert list(block.asu_miller_index) == [(0, 0, 1),
    (1, 0, 0), (1, 0, 0), (2, 0, 0), (2, 0, 0), (2, 2, 2)]

  assert list(block.Ih_values) == [3.0, 3.0, 3.0, 4.0, 4.0, 4.0]
  assert list(block.variances) == [1.0] * 6
  assert list(block.intensities) == [3.0, 1.0, 5.0, 2.0, 6.0, 4.0]

  assert list(block.nonzero_weights) == [2, 0, 4, 1, 5, 3]

  assert block.h_index_matrix[0, 0] == 1
  assert block.h_index_matrix[1, 1] == 1
  assert block.h_index_matrix[2, 1] == 1
  assert block.h_index_matrix[3, 2] == 1
  assert block.h_index_matrix[4, 2] == 1
  assert block.h_index_matrix[5, 3] == 1
  assert block.h_index_matrix.non_zeroes == 6
  assert block.h_index_matrix.n_cols == 4
  assert block.h_index_matrix.n_rows == 6

  # Test for second method to initialise without specifying weights - weights
  # should be set to inverse variances if no weights are given.
  large_reflection_table['variance'] = flex.double(6, 2.0)
  Ih_table = IhTable([(large_reflection_table, None)], test_sg,
    n_blocks=1)
  block = Ih_table.blocked_data_list[0]
  expected_weights = 1.0/flex.double(6, 2.0)
  assert list(block.weights) == list(expected_weights)
  large_reflection_table['variance'] = flex.double(6, 1.0)

  # Test that one can apply an error model, with params that reset w to 1/var
  Ih_table.update_error_model(mock_error_model_2())
  assert list(block.weights) == [1.0, 1.0/2.0, 1.0/3.0, 1.0, 1.0/2.0, 1.0/3.0]

  # Test for functionality of having input selection array.
  sel = flex.bool([True, True, True, False, True, True])
  Ih_table = IhTable([(large_reflection_table, sel)], test_sg,
    weighting_scheme='unity')
  block = Ih_table.blocked_data_list[0]
  assert list(block.intensities) == [3.0, 1.0, 5.0, 2.0, 6.0]
  assert list(block.nonzero_weights) == [2, 0, 4, 1, 5]

  # Expect same behaviour if some are flagged as outliers
  large_reflection_table.set_flags(flex.bool([False, False, False, True, False,
    False]), large_reflection_table.flags.bad_for_scaling)
  Ih_table = IhTable([(large_reflection_table, None)], test_sg,
    weighting_scheme='unity')
  large_reflection_table.set_flags(flex.bool([False, False, False, False, False,
    False]), large_reflection_table.flags.bad_for_scaling)
  block = Ih_table.blocked_data_list[0]
  assert list(block.intensities) == [3.0, 1.0, 5.0, 2.0, 6.0]
  assert list(block.nonzero_weights) == [2, 0, 4, 1, 5]

  # test getting unique group
  table = Ih_table.get_unique_group(flex.miller_index([(1, 0, 0)]))
  assert table.size() == 2
  table = Ih_table.get_unique_group(flex.miller_index([(5, 0, 0)]))
  assert table is None

  # Try to set derivatives and inverse scales
  derivs = sparse.matrix(2, 3)
  derivs[1, 0] = 2.0
  Ih_table.set_derivatives(derivs, 0)
  assert block.derivatives is derivs

  inv_scales = flex.double([1.0, 2.0, 3.0, 4.0, 5.0])
  Ih_table.set_inverse_scale_factors(inv_scales, 0)
  assert list(block.inverse_scale_factors) == list(inv_scales)

  Ih_table.calc_Ih()
  assert list(block.Ih_values) == [3.0, 17.0/13.0, 17.0/13.0, 38.0/41.0,
    38.0/41.0]

def test_IhTable_freework(test_sg):
  """Test the splitting of a multiple dataset into a work and free set."""

  Ih_table = IhTable([(generated_refl_for_splitting_1()[0], None),
    (generated_refl_for_splitting_2()[0], None)], space_group=test_sg,
    free_set_percentage=30.0)

  block = Ih_table.blocked_data_list[0]
  free_block = Ih_table.blocked_data_list[1]

  assert list(block.asu_miller_index) == [(1, 0, 0), (1, 0, 0), (2, 0, 0),
    (2, 0, 0), (1, 0, 0), (2, 0, 0)]
  assert list(free_block.asu_miller_index) == [(0, 0, 1), (2, 2, 2),
    (0, 0, 1), (2, 2, 2), (2, 2, 2)]
  assert Ih_table.free_Ih_table is True
  assert len(Ih_table.blocked_data_list) == 2
  block = Ih_table.blocked_data_list[0]
  free_block = Ih_table.blocked_data_list[1]

  assert block.h_index_matrix[0, 0] == 1
  assert block.h_index_matrix[1, 0] == 1
  assert block.h_index_matrix[2, 1] == 1
  assert block.h_index_matrix[3, 1] == 1
  assert block.h_index_matrix[4, 0] == 1
  assert block.h_index_matrix[5, 1] == 1
  assert block.h_index_matrix.non_zeroes == 6
  assert block.h_index_matrix.n_cols == 2
  assert block.h_index_matrix.n_rows == 6

  assert free_block.h_index_matrix[0, 0] == 1
  assert free_block.h_index_matrix[1, 1] == 1
  assert free_block.h_index_matrix[2, 0] == 1
  assert free_block.h_index_matrix[3, 1] == 1
  assert free_block.h_index_matrix[4, 1] == 1
  assert free_block.h_index_matrix.non_zeroes == 5
  assert free_block.h_index_matrix.n_cols == 2
  assert free_block.h_index_matrix.n_rows == 5

  assert list(Ih_table.blocked_selection_list[0][0]) == [0, 4, 1, 5]
  assert list(Ih_table.blocked_selection_list[0][1]) == [2, 3]
  assert list(Ih_table.blocked_selection_list[1][0]) == [4, 1]
  assert list(Ih_table.blocked_selection_list[1][1]) == [2, 0, 3]

def test_IhTable_split_into_blocks(large_reflection_table,
  small_reflection_table, test_sg):
  """Test that the Ih_table datastructure correctly organises the data
  from two reflection tables into two IhTableBlocks."""

  Ih_table = IhTable([(large_reflection_table, None), (small_reflection_table, None)],
    test_sg, n_blocks=2)

  assert Ih_table.size == 11
  assert Ih_table.n_datasets == 2
  block_list = Ih_table.blocked_data_list

  assert list(block_list[0].miller_index) == [(0, 0, 1), (0, 0, -2),
    (0, 0, 2), (0, 2, 0)]
  assert list(block_list[1].miller_index) == [(0, 4, 0), (1, 0, 0),
    (-1, 0, 0), (0, 4, 0), (0, 4, 0), (1, 0, 0), (10, 0, 0)]

  expected_h_idx_matrix = sparse.matrix(4, 3)
  expected_h_idx_matrix[0, 0] = 1
  expected_h_idx_matrix[1, 1] = 1
  expected_h_idx_matrix[2, 1] = 1
  expected_h_idx_matrix[3, 2] = 1
  assert block_list[0].h_index_matrix == expected_h_idx_matrix
  expected_h_idx_matrix = sparse.matrix(7, 3)
  expected_h_idx_matrix[0, 0] = 1
  expected_h_idx_matrix[1, 1] = 1
  expected_h_idx_matrix[2, 1] = 1
  expected_h_idx_matrix[3, 0] = 1
  expected_h_idx_matrix[4, 0] = 1
  expected_h_idx_matrix[5, 1] = 1
  expected_h_idx_matrix[6, 2] = 1
  assert block_list[1].h_index_matrix == expected_h_idx_matrix

  assert list(block_list[0].Ih_table['dataset_id']) == [0, 0, 0, 0]
  assert list(block_list[1].Ih_table['dataset_id']) == [0, 0, 0, 1, 1, 1, 1]

  block_selection_list = Ih_table.blocked_selection_list

  assert list(block_selection_list[0][0]) == [1, 5, 6, 3]
  assert list(block_selection_list[0][1]) == [4, 0, 2]
  assert list(block_selection_list[1][0]) == []
  assert list(block_selection_list[1][1]) == [1, 3, 0, 2]

def test_set_Ih_values_to_target(test_sg):
  """Test the setting of Ih values for targeted scaling."""
  """Generate input for testing joint_Ih_table."""
  Ih_table = IhTable([(generated_refl_for_splitting_1()[0], None),
    (generated_refl_for_splitting_2()[0], None)], test_sg, n_blocks=2)

  # First check that values are set up correctly.
  block_list = Ih_table.blocked_data_list
  assert list(block_list[0].Ih_values) == [6.0, 17.0/3.0, 17.0/3.0, 6.0, 17.0/3.0]
  assert list(block_list[1].Ih_values) == [16.0/3.0, 16.0/3.0, 7.0, 16.0/3.0, 7.0, 7.0]

  target = IhTable([(generated_refl_for_splitting_1()[0], None),
    (generated_refl_for_splitting_2()[0], None)], test_sg, n_blocks=1)
  #set some values in the target
  target.blocked_data_list[0].Ih_table['Ih_values'] = flex.double([0.1, 0.2,
    0.2, 0.3, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.4])

  for i, block in enumerate(block_list):
    block.set_Ih_values_to_target(target.blocked_data_list[0])
    sel = block.Ih_values != 0.0
    block = block.select(sel)
    Ih_table.apply_selection_to_selection_list(i, sel)
  assert list(block_list[0].Ih_values) == [0.1, 0.2, 0.2, 0.1, 0.2]
  assert list(block_list[1].Ih_values) == [0.3, 0.3, 0.4, 0.3, 0.4, 0.4]
