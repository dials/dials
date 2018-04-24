'''
Tests for outlier rejection.
'''
import pytest
from dials.array_family import flex
from cctbx.sgtbx import space_group
from dials.algorithms.scaling.outlier_rejection import \
  NormDevOutlierRejection, SimpleNormDevOutlierRejection

@pytest.fixture(scope='module')
def test_space_group():
  """Create a space group object."""
  return space_group("C 2y")

@pytest.fixture
def outlier_reflection_table():
  """Create a reflection table with outliers."""
  rt = generate_outlier_table()
  return rt
  
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

def test_standard_outlier_rejection(outlier_reflection_table, test_space_group):
  """Test the outlier rejection algorithm, that the outlier flags are set
  as expected."""
  zmax = 6.0
  refl = NormDevOutlierRejection(outlier_reflection_table,
    test_space_group, zmax).return_reflection_table()
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == [False, False,
    False, False, True, False, False, True, True, True, False, False]

def test_simple_outlier_rejection(outlier_reflection_table, test_space_group):
  """Test the outlier rejection algorithm, that the outlier flags are set
  as expected."""
  zmax = 6.0
  refl = SimpleNormDevOutlierRejection(outlier_reflection_table,
    test_space_group, zmax).return_reflection_table()
  assert list(refl.get_flags(refl.flags.outlier_in_scaling)) == [False, False,
    False, False, True, True, True, True, True, True, False, False]
