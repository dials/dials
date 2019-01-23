"""
Tests for the model_factory module.
"""
import pytest
from mock import Mock
from dials.array_family import flex
from dials.util.options import OptionParser
from dials.algorithms.scaling.model.model import KBScalingModel,\
  PhysicalScalingModel, ArrayScalingModel
from dials.algorithms.scaling.model.scaling_model_factory import \
  KBSMFactory, PhysicalSMFactory, ArraySMFactory, calc_n_param_from_bins
from dials.algorithms.scaling.model.model import initialise_smooth_input
from libtbx import phil

@pytest.fixture
def mock_exp():
  """Return a mock Experiments object."""
  exp = Mock()
  exp.scan.get_oscillation.return_value = [0.0, 1.0]
  exp.scan.get_oscillation_range.return_value = [0, 90]
  return exp

@pytest.fixture(scope='module')
def default_params():
  """Return the default parsed params phil scope."""
  return generated_param()

@pytest.fixture
def test_reflections():
  """Return a reflection table"""
  return generated_refl()

@pytest.fixture()
def mock_physical_params():
  """Return a mock params object for a physical model."""
  params = Mock()
  params.parameterisation.scale_term = True
  params.parameterisation.scale_interval = 10.0
  params.parameterisation.decay_term = True
  params.parameterisation.decay_interval = 15.0
  params.parameterisation.absorption_term = True
  params.parameterisation.lmax = 4
  params.parameterisation.decay_restraint = 1e-1
  return params

def generated_refl():
  """Create a reflection table."""
  rt = flex.reflection_table()
  rt['d'] = flex.double([1.0, 1.0, 1.0, 1.0])
  rt['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  rt.set_flags(flex.bool([False, False, False, False]),
    rt.flags.user_excluded_in_scaling)
  return rt

def generated_param():
  """Generate the default scaling parameters object."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  parameters.parameterisation.modulation_term = True
  return parameters

def test_ScalingModelfactories(default_params, mock_exp, test_reflections):
  """Test the factory creation of the three standard scaling models with the
  default params."""

  KBmodel = KBSMFactory.create(default_params, [], [])
  assert isinstance(KBmodel, KBScalingModel)

  physicalmodel = PhysicalSMFactory.create(default_params, mock_exp, test_reflections)
  assert isinstance(physicalmodel, PhysicalScalingModel)

  arraymodel = ArraySMFactory.create(default_params, mock_exp, test_reflections)
  assert isinstance(arraymodel, ArrayScalingModel)

  # Add more rigorous tests to checl that the model has been set up correctly.?
  # Might be best to refactor scaling model factories first.


def test_PhysicalSMFactory(mock_physical_params, mock_exp, test_reflections):
  """Test that it passes the correct dict to physical model."""
  physicalmodel = PhysicalSMFactory.create(mock_physical_params, mock_exp,
    test_reflections)
  assert isinstance(physicalmodel, PhysicalScalingModel)
  assert physicalmodel.configdict['lmax'] == (
    mock_physical_params.parameterisation.lmax)
  assert physicalmodel.components['absorption'].n_params == 24
  assert list(physicalmodel.components['absorption'].parameters) == [0.0] * 24

def test_model_factory_utilities(mock_exp):
  """Test the utility functions in the scaling_model_factory module."""

  # Test calc_n_param_from_bins(value_min, value_max, n_bins)
  assert calc_n_param_from_bins(0.0, 1.0, 1) == (2, 1.0)
  assert calc_n_param_from_bins(0.0, 2.0, 2) == (3, 1.0)
  assert calc_n_param_from_bins(0.0, 3.0, 3) == (5, 1.0)
  assert calc_n_param_from_bins(0.0, 10.0, 10) == (12, 1.0)
  assert calc_n_param_from_bins(0.0, 10.0, 5) == (7, 2.0)
  with pytest.raises(AssertionError):
    (_, _) = calc_n_param_from_bins(0.0, 1.0, 0)
  with pytest.raises(AssertionError):
    (_, _) = calc_n_param_from_bins(0.0, 1.0, 0.5)

  # Test initialise_smooth_input(osc_range, one_osc_width, interval)
  # This is initialised with the oscillation range, width of one osc and
  # rotation interval in degress, returning
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 1.0, 1.0)
  assert (n_param, norm_fac, rot_int) == (12, 0.999, 1.0)
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 1.0, 12)
  assert (n_param, norm_fac, rot_int) == (2, 0.0999, 10.0)
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 1.0, 10)
  assert (n_param, norm_fac, rot_int) == (2, 0.0999, 10.0)
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 1.0, 9.99)
  assert (n_param, norm_fac, rot_int) == (3, 0.1998, 5.0)
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 1.0, 5.0)
  assert (n_param, norm_fac, rot_int) == (3, 0.1998, 5.0)
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 1.0, 4.99)
  assert (n_param, norm_fac, rot_int) == (5, 0.999*3.0/10.0, 10.0/3.0)
  n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 2.0, 4.99)
  assert (n_param, norm_fac, rot_int) == (5, 0.999*6.0/10.0, 10.0/3.0)
