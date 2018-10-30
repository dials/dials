"""
Tests for the post_scaling_analysis module.
"""
import pytest
from mock import Mock
from dxtbx.model import Crystal
from dials.array_family import flex
from dials.algorithms.scaling.post_scaling_analysis import \
  exclude_on_batch_rmerge, exclude_on_image_scale
from dials.algorithms.scaling.model.model import PhysicalScalingModel, KBScalingModel


@pytest.fixture
def mock_exp():
  """Mock experiment object for tests."""
  exp = Mock()
  parameters_dict = {'scale': {'parameters' : flex.double([1.0, 1.0, 0.6, 0.3, 0.2]),
    'parameter_esds' : flex.double([0.1, 0.1, 0.1, 0.1, 0.1])}}
  configdict = {'corrections' : ['scale'], 's_norm_fac' : 1.0}
  exp.scaling_model = PhysicalScalingModel(parameters_dict, configdict)
  exp.scan.get_image_range.return_value = [1, 3]
  exp.scan.get_oscillation.return_value = [0, 1]
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  exp.crystal = Crystal.from_dict(exp_dict)
  return exp

@pytest.fixture
def mock_2_KB_exps():
  """Mock experiments with scaling models for tests."""
  exps = []
  exp = Mock()
  parameters_dict_1 = {'scale': {'parameters' : flex.double([0.2]),
    'parameter_esds' : flex.double([0.1])}}
  configdict = {'corrections' : ['scale']}
  exp.scaling_model = KBScalingModel(parameters_dict_1, configdict)
  exps.append(exp)
  exp = Mock()
  parameters_dict_2 = {'scale': {'parameters' : flex.double([1.0]),
    'parameter_esds' : flex.double([0.1])}}
  exp.scaling_model = KBScalingModel(parameters_dict_2, configdict)
  exps.append(exp)
  return exps

def generate_test_refl():
  """Make a test reflection table."""
  refl = flex.reflection_table()
  refl['intensity'] = flex.double([1.0, 0.2, 0.3])
  refl['xyzobs.px.value'] = flex.vec3_double([(0, 0, 0.5), (0, 0, 1.5), (0, 0, 2.5)])
  refl.set_flags(flex.bool(refl.size(), False), refl.flags.user_excluded_in_scaling)
  return refl

def generate_bad_rmerge_refl():
  """Make a test reflection table."""
  refl = flex.reflection_table()
  refl['intensity'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0])
  refl['variance'] = flex.double([1.0] * 9 + [100])
  refl['miller_index'] = flex.miller_index([(1, 0, 0)] * 10)
  refl['xyzobs.px.value'] = flex.vec3_double([(0, 0, 0.5)] * 9 + [(0, 0, 2.5)])
  refl.set_flags(flex.bool(refl.size(), False), refl.flags.user_excluded_in_scaling)
  return refl

@pytest.fixture
def test_2_refls():
  """Make a list of two reflection tables."""
  refls = []
  refl = flex.reflection_table()
  refl['intensity'] = flex.double([1.0, 0.2, 0.3])
  refl.set_flags(flex.bool(refl.size(), False), refl.flags.user_excluded_in_scaling)
  refls.append(refl)
  refl = flex.reflection_table()
  refl['intensity'] = flex.double([1.0, 0.2, 0.3])
  refl.set_flags(flex.bool(refl.size(), False), refl.flags.user_excluded_in_scaling)
  refls.append(refl)
  return refls

def test_exclude_on_image_scale(mock_exp, mock_2_KB_exps, test_2_refls):
  """Test excluding certain images based on the inverse scale factor of the image."""
  # Set up an example with physical scaling model.
  assert mock_exp.scaling_model.id_ == 'physical'
  # Need to configure the reflection table and model
  # The scale component is parameterised by 5 values, at normalised rotation of
  # -0.5, 0.5, 1.5, 2.5 3.5 (i.e. for data spanning 0>3 normalised rot)
  test_refl = generate_test_refl()
  mock_params = Mock()
  mock_params.parameterisation.decay_restraint = 0.0
  test_refl = mock_exp.scaling_model.configure_reflection_table(test_refl, mock_exp, mock_params)
  mock_exp.scaling_model.components['scale'].update_reflection_data(test_refl)
  new_refl = exclude_on_image_scale([test_refl], [mock_exp], 0.5)[0]
  # Excluding below 0.5 should only pick out the last reflection.

  assert list(new_refl.get_flags(new_refl.flags.user_excluded_in_scaling)) == [
    False, False, True]

  # Test again with a KB model dataset with no scan object - should exclude
  # all of the first dataset and none of the second.
  new_refls = exclude_on_image_scale(test_2_refls, mock_2_KB_exps, 0.5)
  assert list(new_refls[0].get_flags(
    new_refls[0].flags.user_excluded_in_scaling)) == [True, True, True]
  assert list(new_refls[1].get_flags(
    new_refls[1].flags.user_excluded_in_scaling)) == [False, False, False]

  # Test again on another kind of scaling model - should gracefully do nothing
  exp = Mock()
  exp.scaling_model.id_ = 'New'
  test_refl = generate_test_refl()
  new_refl = exclude_on_image_scale([test_refl], [exp], 0.5)[0]
  assert list(new_refl.get_flags(new_refl.flags.user_excluded_in_scaling)) == [
    False, False, False]

def test_exclude_on_batch_rmerge(mock_exp):
  """Test excluding data based on batch rmerge."""
  # Want a dataset with a bad r-merge for one or two frames -
  # set low to trigger finding in last frame
  refl = generate_bad_rmerge_refl()
  new_refl = exclude_on_batch_rmerge([refl], [mock_exp], 0.5)[0]
  assert list(new_refl.get_flags(new_refl.flags.user_excluded_in_scaling)) == [
    False] * 9 + [True]

  # Test again but for no scan - make a good and a bad frame and a frame with
  # no good data
  mock_exp.scan = False
  refl1 = generate_bad_rmerge_refl()
  refl1['intensity'][9] = 1.0
  refl2 = generate_bad_rmerge_refl()
  refl3 = generate_bad_rmerge_refl()
  refl3.set_flags(flex.bool(refl3.size(), True), refl3.flags.outlier_in_scaling)
  new_refl = exclude_on_batch_rmerge([refl1, refl2, refl3], [mock_exp, mock_exp, mock_exp], 0.5)
  assert list(new_refl[0].get_flags(new_refl[0].flags.user_excluded_in_scaling)) == [
    False] * 10
  assert list(new_refl[1].get_flags(new_refl[1].flags.user_excluded_in_scaling)) == [
    True] * 10
  assert list(new_refl[2].get_flags(new_refl[2].flags.user_excluded_in_scaling)) == [
    False] * 10
