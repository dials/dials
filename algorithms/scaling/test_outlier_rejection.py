'''
Tests for outlier rejection.
'''
import pytest
from libtbx.utils import Sorry
from cctbx.sgtbx import space_group
from dials.array_family import flex
from dials.algorithms.scaling.outlier_rejection import \
  NormDevOutlierRejection, SimpleNormDevOutlierRejection, reject_outliers,\
  TargetedOutlierRejection

@pytest.fixture(scope='module')
def test_sg():
  """Create a space group object."""
  return space_group("C 2y")

@pytest.fixture
def outlier_reflection_table():
  """Create a reflection table with outliers."""
  rt = generate_outlier_table()
  return rt

@pytest.fixture
def outlier_target_table():
  target = flex.reflection_table()
  target['intensity'] = flex.double([500])
  target['variance'] = flex.double([1.0])
  target['inverse_scale_factor'] = flex.double([1.0])
  target['miller_index'] = flex.miller_index([(0, 0, 2)])
  target.set_flags(flex.bool([False]), target.flags.excluded_for_scaling)
  target.set_flags(flex.bool([False]), target.flags.user_excluded_in_scaling)
  return target

def generate_outlier_table():
  """Generate a reflection table for outlier testing."""
  rt = flex.reflection_table()
  rt['intensity'] = flex.double([1.0, 1.0, 1.0, 1.0, 20.0, 1.0, 1.0, 20.0,
    400.0, 500.0, 10.0, 10.0])
  rt['variance'] = flex.double(12, 1.0)
  rt['inverse_scale_factor'] = flex.double(12, 1.0)
  rt['miller_index'] = flex.miller_index([(0, 0, 1), (0, 0, 1), (0, 0, 1),
    (0, 0, 1), (0, 0, 1), (0, 0, 2), (0, 0, 2), (0, 0, 2), (0, 0, 2), (0, 0, 2),
    (0, 0, 23), (0, 0, 3)])
  rt.set_flags(flex.bool([True, False, False, False, False, False, False, False,
    False, False, False, False]),
    rt.flags.excluded_for_scaling)
  rt.set_flags(flex.bool(12, False), rt.flags.user_excluded_in_scaling)
  return rt

expected_standard_output = [False, False, False, False, True, False, False,
  True, True, True, False, False]

expected_simple_output = [False, False, False, False, True, True, True, True,
  True, True, False, False]

expected_target_output = [False, False, False, False, False, True, True, True,
  True, False, False, False]

def test_standard_outlier_rejection(outlier_reflection_table, test_sg):
  """Test the outlier rejection algorithm, that the outlier flags are set
  as expected."""
  zmax = 6.0
  OutlierRej = NormDevOutlierRejection(outlier_reflection_table,
    test_sg, zmax)
  refl = OutlierRej.return_reflection_table()
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == \
    expected_standard_output
  OutlierRej.unset_outlier_flags()
  refl = OutlierRej.return_reflection_table()
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == [False] * 12

def test_targeted_outlier_rejection(outlier_reflection_table,
  outlier_target_table, test_sg):
  """Test the targeted outlier rejection algorithm - only reflections
  that exist in both the target and the reflecton table should be tested
  based on their normalised deviation."""
  zmax = 6.0
  refl = TargetedOutlierRejection(outlier_reflection_table,
    test_sg, zmax, outlier_target_table).return_reflection_table()
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == \
    expected_target_output

def test_simple_outlier_rejection(outlier_reflection_table, test_sg):
  """Test the outlier rejection algorithm, that the outlier flags are set
  as expected."""
  zmax = 6.0
  refl = SimpleNormDevOutlierRejection(outlier_reflection_table,
    test_sg, zmax).return_reflection_table()
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == \
    expected_simple_output

def test_reject_outliers_function(outlier_reflection_table,
  outlier_target_table, test_sg):
  """Test the helper function."""

  refl = reject_outliers(outlier_reflection_table, test_sg, 'standard')
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == \
    expected_standard_output

  refl = reject_outliers(outlier_reflection_table, test_sg, 'simple')
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == \
    expected_simple_output

  refl = reject_outliers(outlier_reflection_table, test_sg, 'target',
    target=outlier_target_table)
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == \
    expected_target_output

  with pytest.raises(Sorry):
    refl = reject_outliers(outlier_reflection_table, test_sg, 'badchoice')
