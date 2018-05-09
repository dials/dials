"""
Tests for the scaler factory classes and helper functions.
"""
from copy import deepcopy
import pytest
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.utils import Sorry
from dxtbx.model import Crystal
from mock import Mock, MagicMock
from dials.algorithms.scaling.scaler_factory import SingleScalerFactory,\
  TargetScalerFactory, MultiScalerFactory, is_scaled, create_scaler
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  MultiScaler, TargetScaler, NullScaler

def generated_refl():
  """Generate a test reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['inverse_scale_factor'] = flex.double(4, 1.0)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (2, 0, 0), (2, 2, 2)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0, 0.0]) #don't change
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, False, False]),
    reflections.flags.integrated)
  reflections.set_flags(flex.bool([False, False, True, True]),
    reflections.flags.bad_for_scaling)
  return reflections

@pytest.fixture
def test_refl():
  """Generate a test reflection table."""
  return generated_refl()

@pytest.fixture
def refl_list():
  """Make a list of three reflection tables."""
  refl_list = [generated_refl()]
  refl_list.append(generated_refl())
  refl_list.append(generated_refl())
  return refl_list

@pytest.fixture
def generated_param():
  """Generate a param phil scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  parameters.scaling_options.free_set_percentage = 50.0
  return parameters

@pytest.fixture
def mock_scaling_component():
  """Mock scaling component to allow creation of a scaling model."""
  component = MagicMock()
  component.n_params = 2
  component.inverse_scales = flex.double([0.9, 1.1])
  return component

@pytest.fixture
def mock_exp(mock_scaling_component):
  """Mock experiments object for initialising a scaler."""

  def side_effect_config_table(*args):
    """Side effect to mock configure reflection table
    call during initialisation."""
    return args[0]

  exp = Mock()
  exp.scaling_model.components = {'scale' : mock_scaling_component}
  exp.scaling_model.consecutive_refinement_order = ['scale']
  exp.scaling_model.is_scaled = False
  exp.scaling_model.configure_reflection_table.side_effect = side_effect_config_table
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  exp.crystal = Crystal.from_dict(exp_dict)
  return exp

@pytest.fixture
def mock_explist_3exp(mock_exp):
  """A mock experimentlist, containing one mock exp instance three times."""
  exp = [mock_exp]
  exp.append(mock_exp)
  exp.append(mock_exp)
  return exp

@pytest.fixture
def mock_scaled_exp():
  """A mock experiments object with scaling_model.is_scaled = True"""
  exp = Mock()
  exp.scaling_model.is_scaled = True
  return exp

@pytest.fixture
def mock_unscaled_exp():
  """A mock experiments object with scaling_model.is_scaled = False"""
  exp = Mock()
  exp.scaling_model.is_scaled = False
  return exp

@pytest.fixture
def mock_experimentlist(mock_scaled_exp, mock_unscaled_exp):
  """A mock experimentlist of mock scaled/unscaled mock exp."""
  explist = [mock_scaled_exp, mock_scaled_exp,
    mock_unscaled_exp, mock_scaled_exp, mock_unscaled_exp]
  return explist

def test_SingleScalerFactory(generated_param, mock_exp, test_refl):
  """Test the single scaler factory."""

  #Test default, (id = 0, no split into free set)
  ss = SingleScalerFactory.create(generated_param, mock_exp, test_refl)
  assert isinstance(ss, SingleScalerBase)
  assert set(ss.reflection_table['id']) == set([0])
  assert list(ss.Ih_table.free_set_sel) == []

  # Now set the free set option and give an id
  generated_param.scaling_options.use_free_set = True
  ss = SingleScalerFactory.create(generated_param, mock_exp, test_refl,
    scaled_id=1)
  assert set(ss.reflection_table['id']) == set([1])
  assert list(ss.Ih_table.free_set_sel) != []

  # Nowflag that the single scaler is being created for a multiscaler -
  # therefore the Ih table should not be split into a free set.
  ss = SingleScalerFactory.create(generated_param, mock_exp, test_refl,
    scaled_id=1, for_multi=True)
  assert list(ss.Ih_table.free_set_sel) == []

def test_TargetScalerFactory(generated_param, mock_explist_3exp, refl_list):
  """Test the target scaler factory."""

  # Test standard initialisation.
  assert generated_param.scaling_options.use_free_set is False #just to check
  target = TargetScalerFactory.create(generated_param, mock_explist_3exp,
    refl_list, is_scaled_list=[True, True, False])
  assert isinstance(target, TargetScaler)
  assert len(target.single_scalers) == 2
  assert len(target.unscaled_scalers) == 1
  assert set(target.single_scalers[0].reflection_table['id']) == set([0])
  assert set(target.single_scalers[1].reflection_table['id']) == set([1])
  assert set(target.unscaled_scalers[0].reflection_table['id']) == set([2])
  assert list(target.single_scalers[0].Ih_table.free_set_sel) == []
  assert list(target.single_scalers[1].Ih_table.free_set_sel) == []

  # Test initialisation with free set option - still should not cause free
  # set splitting in single scalers.
  generated_param.scaling_options.use_free_set = True
  target = TargetScalerFactory.create(generated_param, mock_explist_3exp,
    refl_list, is_scaled_list=[True, True, False])
  assert list(target.single_scalers[0].Ih_table.free_set_sel) == []
  assert list(target.single_scalers[1].Ih_table.free_set_sel) == []

  # Test for correct initialisation hen scaling against a target model.
  generated_param.scaling_options.target_model = True
  target = TargetScalerFactory.create(generated_param, mock_explist_3exp,
    refl_list, is_scaled_list=[True, True, False])
  assert isinstance(target.single_scalers[0], NullScaler)
  assert isinstance(target.single_scalers[1], NullScaler)

  # Now test converting targetscaler to multiscaler
  multiscaler = MultiScalerFactory.create_from_targetscaler(target)
  assert isinstance(multiscaler, MultiScaler)
  assert len(multiscaler.single_scalers) == 3

def test_MultiScalerFactory(generated_param, mock_explist_3exp, refl_list):
  """Test the MultiScalerFactory."""

  multiscaler = MultiScalerFactory.create(generated_param, mock_explist_3exp,
    refl_list)
  assert isinstance(multiscaler, MultiScaler)
  assert len(multiscaler.single_scalers) == 3
  for i in range(3):
    assert set(multiscaler.single_scalers[i].reflection_table['id']) == set([i])

  generated_param.scaling_options.use_free_set = True
  multiscaler = MultiScalerFactory.create(generated_param, mock_explist_3exp,
    refl_list)
  assert list(multiscaler.single_scalers[0].Ih_table.free_set_sel) == []

def test_scaler_factory_helper_functions(mock_experimentlist, mock_exp,
  mock_explist_3exp, generated_param, test_refl, refl_list):
  """Test the helper functions."""

  #Test is_scaled function
  scaled_list = is_scaled(mock_experimentlist)
  assert scaled_list == [True, True, False, True, False]

  #Test create_scaler
  #Test case for single refl and exp
  scaler = create_scaler(generated_param, [mock_exp], [test_refl])
  assert isinstance(scaler, SingleScalerBase)

  #If none or allscaled
  scaler = create_scaler(generated_param, mock_explist_3exp, refl_list)
  assert isinstance(scaler, MultiScaler)

  mock_explist_3exp[0].scaling_model.is_scaled = False
  # ^ changes all in list as same instance of exp.
  scaler = create_scaler(generated_param, mock_explist_3exp, refl_list)
  assert isinstance(scaler, MultiScaler)

  #If only some scaled
  mock_explist_3exp[1] = deepcopy(mock_explist_3exp[0])
  mock_explist_3exp[1].scaling_model.is_scaled = True
  scaler = create_scaler(generated_param, mock_explist_3exp, refl_list)
  assert isinstance(scaler, TargetScaler)

  #If no reflections passed in.
  with pytest.raises(Sorry):
    scaler = create_scaler(generated_param, mock_explist_3exp, [])
