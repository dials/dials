'''
This code tests for Ih_table and joint_Ih_table data structures.
This also provides a test for the scaler, which must be successfully
initialised in order to provide input for the Ih_table.
'''
from copy import deepcopy
import pytest
from dials.algorithms.scaling.Ih_table import IhTable
from dials.array_family import flex
from cctbx.sgtbx import space_group
from scitbx import sparse

@pytest.fixture()
def reflection_table_for_block():
  return generated_refl_for_splitting_1()

@pytest.fixture()
def large_reflection_table():
  """Create a reflection table."""
  return generate_refl_1()

def generated_refl_for_splitting_1():
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

def generated_refl_for_splitting_2():
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

@pytest.fixture
def joint_test_input(large_reflection_table, small_reflection_table, test_sg):
  """Generate input for testing joint_Ih_table."""
  Ih_table_1 = IhTable([(large_reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')
  Ih_table_2 = IhTable([(small_reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')
  return Ih_table_1, Ih_table_2

def test_Ih_table(reflection_table_for_block, test_sg):
  """Test for Ih_table datastructure. Upon initialisation, Ih_table should set
  unity scale factors and calculate Ih_values. It should also create the
  a h_index_matrix."""
  large_reflection_table = reflection_table_for_block[0]
  large_reflection_table['variance'] = flex.double(6, 2.0)
  # Test initialisation fails without correct columns
  for key in ['intensity', 'variance', 'inverse_scale_factor', 'miller_index']:
    col_copy = deepcopy(large_reflection_table[key])
    del large_reflection_table[key]
    with pytest.raises(AssertionError):
      Ih_table = IhTable([(large_reflection_table, None)], test_sg,
        n_blocks=1, weighting_scheme='unity')
    large_reflection_table[key] = col_copy

  weights = flex.double(6, 1.0)
  Ih_table = IhTable([(large_reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')

  assert Ih_table.id_ == "IhTable"

  # Tests calc_Ih, assign_h_matrices, interface
  block = Ih_table.blocked_data_list[0]
  assert block.size == 6
  assert (block.inverse_scale_factors == 1.0).count(True) == block.size
  assert (block.weights == 1.0).count(True) == block.size
  assert list(block.asu_miller_index) == list(flex.miller_index([(0, 0, 1),
    (1, 0, 0), (1, 0, 0), (2, 0, 0), (2, 0, 0), (2, 2, 2)]))

  assert list(block.Ih_values) == list(flex.double(
    [3.0, 3.0, 3.0, 4.0, 4.0, 4.0]))
  assert list(block.variances) == [2.0] * 6
  assert list(block.intensities) == list(flex.double(
    [3.0, 1.0, 5.0, 2.0, 6.0, 4.0]))

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

  # Test calc_nh function.
  block.calc_nh()
  assert list(block.n_h) == [1.0, 2.0, 2.0, 2.0, 2.0, 1.0]

  # Test selection function
  sel = flex.bool([True, True, False, False, False, False])
  block.select(sel)
  assert block.size == 2
  assert (block.inverse_scale_factors == 1.0).count(True) == block.size
  assert (block.weights == 1.0).count(True) == block.size
  assert list(block.asu_miller_index) == list(flex.miller_index(
    [(0, 0, 1), (1, 0, 0)]))
  assert list(block.Ih_values) == list(flex.double([3.0, 3.0]))
  assert list(block.variances) == list(flex.double([2.0, 2.0]))
  assert list(block.intensities) == list(flex.double([3.0, 1.0]))
  assert block.h_index_matrix[0, 0] == 1
  assert block.h_index_matrix[1, 1] == 1
  assert block.h_index_matrix.non_zeroes == 2
  assert block.h_index_matrix.n_cols == 2
  assert block.h_index_matrix.n_rows == 2
  # Test that non-zero weights was correctly updated
  assert list(block.nonzero_weights) == [2, 0]

  # Test for second method to initialise without specifying weights - weights
  # should be set to inverse variances if no weights are given.
  Ih_table = IhTable([(large_reflection_table, None)], test_sg,
    n_blocks=1)
  block = Ih_table.blocked_data_list[0]
  expected_weights = 1.0/flex.double(6, 2.0)
  assert list(block.weights) == list(expected_weights)
  # Test that one can set the weights to unity
  block.weights = weights
  assert (block.weights == 1.0).count(True) == block.size

  # Test that one can apply an error model, with params that reset w to 1/var
  #Ih_table.update_error_model([1.0, 0.0])
  #assert approx_equal(list(block.weights), list(expected_weights))

  # note - is the below feature actually needed anymore?
  '''# Test for functionality of having preset Ih_values, set to a tenth of what
  # they would be. Test that these are set after Ih table is initialised.
  large_reflection_table['Ih_values'] = flex.double(
    [0.3, 0.4, 0.3, 0.4, 0.3, 0.4])
  #Ih_table = SingleIhTable(reflection_table, test_sg, weighting_scheme='unity', split_blocks=None)
  Ih_table = IhTable([(large_reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')
  block = Ih_table.blocked_data_list[0]
  assert list(block.Ih_values) == list(flex.double(
    [0.3, 0.3, 0.3, 0.4, 0.4, 0.4]))'''

  # Test for size checking when setting new weights, inverse scales
  with pytest.raises(AssertionError):
    block.inverse_scale_factors = flex.double([1.0])
  with pytest.raises(AssertionError):
    block.weights = flex.double([1.0])

  # test getting unique group
  table = Ih_table.get_unique_group(flex.miller_index([(1, 0, 0)]))
  assert table.size() == 2
  table = Ih_table.get_unique_group(flex.miller_index([(5, 0, 0)]))
  assert table is None

@pytest.mark.skip(reason='Free set selection not yet implemented.')
def test_Ih_table_freework(large_reflection_table, test_sg):
  """Test the splitting of a single dataset into a work and free set."""
  reflection_table = large_reflection_table
  Ih_table = IhTable([(reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')

  old_Ih_values = deepcopy(Ih_table.Ih_values)

  Ih_table.split_into_free_work(30.0)

  #Test that the Ih_table was split properly.
  assert Ih_table.free_Ih_table
  assert list(Ih_table.free_set_sel) == [False, True, False, False, True,
    False, False]
  assert Ih_table.free_Ih_table.h_index_matrix.non_zeroes == 2
  assert Ih_table.free_Ih_table.h_index_matrix.n_cols == 2
  assert Ih_table.free_Ih_table.h_index_matrix[0, 0] == 1
  assert Ih_table.free_Ih_table.h_index_matrix[1, 1] == 1

  assert Ih_table.h_index_matrix.non_zeroes == 5
  assert Ih_table.h_index_matrix.n_cols == 3
  assert Ih_table.h_index_matrix[0, 2] == 1
  assert Ih_table.h_index_matrix[1, 2] == 1
  assert Ih_table.h_index_matrix[2, 1] == 1
  assert Ih_table.h_index_matrix[3, 0] == 1
  assert Ih_table.h_index_matrix[4, 0] == 1

  # Now test that calc_Ih correctly calculates for the two sets.
  Ih_table.calc_Ih()
  assert list(Ih_table.Ih_values) == list(
    old_Ih_values.select(~Ih_table.free_set_sel))
  assert list(Ih_table.free_Ih_table.Ih_values) == list(
    old_Ih_values.select(Ih_table.free_set_sel))

  # Now test setting of new inverse scales in the case of a free set
  new_scales = flex.double([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
  Ih_table.inverse_scale_factors = [new_scales]
  assert Ih_table.inverse_scale_factors == new_scales.select(
    ~Ih_table.free_set_sel)
  assert Ih_table.free_Ih_table.inverse_scale_factors == new_scales.select(
    Ih_table.free_set_sel)

@pytest.mark.skip(reason='Free set selection not yet implemented.')
def test_joint_Ih_table_freework(large_reflection_table, small_reflection_table, test_sg):
  """Test the splitting of a multiple dataset into a work and free set."""
  #Ih_table_1, Ih_table_2 = joint_test_input
  Ih_table = IhTable([(large_reflection_table[0], None),
    (small_reflection_table[0], None)], space_group=test_sg)
  old_inverse_scales = Ih_table.inverse_scale_factors
  old_Ih_vals = Ih_table.Ih_values
  Ih_table.split_into_free_work(20.0)
  # As the split method is contained in the base, shouldn't need to check all
  # elements, but do need to check the calc_Ih method which is a bit different.
  assert Ih_table.free_set_sel
  assert Ih_table.h_index_matrix.non_zeroes == 9
  assert Ih_table.h_index_matrix.n_cols == 4
  assert Ih_table.h_index_matrix.n_rows == 9
  assert Ih_table.free_Ih_table.h_index_matrix.non_zeroes == 2
  assert Ih_table.free_Ih_table.h_index_matrix.n_cols == 2
  assert Ih_table.free_Ih_table.h_index_matrix.n_rows == 2
  Ih_table.calc_Ih()
  assert list(Ih_table.inverse_scale_factors) == list(
    old_inverse_scales.select(~Ih_table.free_set_sel))
  assert list(Ih_table.Ih_values) == list(
    old_Ih_vals.select(~Ih_table.free_set_sel))
  assert list(Ih_table.free_Ih_table.inverse_scale_factors) == list(
    old_inverse_scales.select(Ih_table.free_set_sel))
  assert list(Ih_table.free_Ih_table.Ih_values) == list(
    old_Ih_vals.select(Ih_table.free_set_sel))

def test_Ih_table_nonzero_weights(large_reflection_table, test_sg):
  """Test for 'nonzero_Weights' attribute and how this changes during selection.
  The purpose of this is to indicate the relationship of the Ih_table data to
  the original input reflection table."""
  reflection_table = large_reflection_table
  Ih_table = IhTable([(reflection_table, None)], test_sg,
    n_blocks=1, weighting_scheme='unity')
  block = Ih_table.blocked_data_list[0]
  assert list(block.nonzero_weights) == [1, 5, 6, 3, 4, 0, 2]
  block.select(flex.bool([True, True, True, False, False, False,
    False]))
  assert block.size == 3
  assert list(block.nonzero_weights) == [1, 5, 6]
  reflection_table.set_flags(flex.bool([True, False, False, False, False,
    False, False]), reflection_table.flags.bad_for_scaling)
  Ih_table = IhTable([(reflection_table, None)], test_sg,
    n_blocks=1)
  block = Ih_table.blocked_data_list[0]
  assert block.size == 6
  assert list(block.nonzero_weights) == [1, 5, 6, 3, 4, 2]
  block.select(flex.bool([True, True, True, False, False, False]))
  assert list(block.nonzero_weights) == [1, 5, 6]

def test_iterative_weighting(large_reflection_table, test_sg):
  """Test the setting of iterative weights."""
  rt = large_reflection_table
  Ih_table = IhTable([(rt, None)], test_sg,
    n_blocks=1, weighting_scheme='GM')
  #Ih_table = SingleIhTable(rt, test_sg, weighting_scheme='GM', split_blocks=None)

  # After construction, weights should initially be set to inverse variances,
  # to allow the first calculation of Ih (by least-squares approach).
  assert list(Ih_table.blocked_data_list[0].weights) == [1.0/100.0, 1.0/50.0,
    1.0/50.0, 1.0/60.0, 1.0/30.0, 1.0/90.0, 1.0/90.0]

  #  list(1.0/rt['variance'])

  # Now update weights
  Ih_table.update_weights()
  block = Ih_table.blocked_data_list[0]
  gIh = block.inverse_scale_factors * block.Ih_values
  t = (block.intensities - gIh) / gIh
  assert list(block.weights) == list(1.0/(1.0 + t**2)**2)


def test_set_Ih_values_to_target(test_sg):
  """Test the setting of Ih values for targeted scaling."""
  """Generate input for testing joint_Ih_table."""
  Ih_table = IhTable([(generated_refl_for_splitting_1()[0], None),
    (generated_refl_for_splitting_2()[0], None)], test_sg, n_blocks=2)
  #Ih_table_1 = SingleIhTable(large_reflection_table, test_sg,
  #  weighting_scheme='unity', split_blocks=None)
  #Ih_table_2 = SingleIhTable(small_reflection_table, test_sg,
  #  weighting_scheme='unity', split_blocks=None)
  #target_Ih_table, Ih_table_2 = joint_test_input

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


def test_new_joint_Ih_table(large_reflection_table, small_reflection_table, test_sg):
  """Test that the joint_Ih_table datastructure correctly combined the data
  from two reflection tables."""

  Ih_table = IhTable([(large_reflection_table, None), (small_reflection_table, None)],
    test_sg, n_blocks=2)

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

  # Test setting of error models and updating weights.
  '''Ih_table_1.update_error_model([0.5, 0.0])
  Ih_table_2.update_error_model([0.5, 0.0])
  Ih_table.update_weights_from_error_models()
  expected_weights = 4.0/flex.double(
    [90.0, 100.0, 90.0, 60.0, 30.0, 50.0, 50.0, 60.0, 30.0, 10.0, 30.0])
  assert approx_equal(list(Ih_table.weights), list(expected_weights))'''

def test_split_into_blocks(large_reflection_table, test_sg):
  """Test the setting of iterative weights."""
  rt = large_reflection_table
  rt2 = deepcopy(large_reflection_table)
  rt.extend(rt2)
  _ = IhTable([(rt, None)], test_sg, n_blocks=2)
