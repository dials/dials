'''
This code tests the data managers and active parameter managers.
'''
import pytest
import mock
from mock import Mock, MagicMock, call
from scitbx import sparse
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from cctbx.sgtbx import space_group
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.scaling_library import create_scaling_model
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  calc_sf_variances, ScalerBase, MultiScalerBase, MultiScaler, TargetScaler
from dials.algorithms.scaling.Ih_table import JointIhTable, SingleIhTable

@pytest.fixture
def test_reflections():
  """Make a test reflection table."""
  return generated_refl()

@pytest.fixture
def test_reflections_Ihtable():
  """Make a test reflection table."""
  return generated_refl_Ihtable()

@pytest.fixture
def test_experiments():
  """Make a test experiments list"""
  return generated_exp()

@pytest.fixture
def test_params():
  """Make a test param phil scope."""
  return generated_param()

@pytest.fixture
def test_sg():
  """Make a test param phil scope."""
  return space_group("C 2y")

def generated_refl_Ihtable():
  """A test reflection tale for instantiating a SingleIhTable"""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['variance'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['inverse_scale_factor'] = flex.double(4, 1.0)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (2, 0, 0), (2, 2, 2)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0, 0.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, False, False]),
    reflections.flags.integrated)
  reflections.set_flags(flex.bool([False, False, True, True]),
    reflections.flags.bad_for_scaling)
  return [reflections]

def generated_refl():
  """Generate a reflection table."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['intensity.sum.value'] = flex.double([10.0, 100.0, 1000.0, 10.0])
  reflections['intensity.sum.variance'] = flex.double([10.0, 100.0, 1000.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (2, 0, 0), (2, 2, 2)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0, 0.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, False, False]),
    reflections.flags.integrated)
  reflections.set_flags(flex.bool([False, False, True, True]),
    reflections.flags.bad_for_scaling)
  return [reflections]

def generated_exp():
  """Generate an experiment list with two experiments."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

def generated_param():
  """Generate a param phil scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  parameters.scaling_options.space_group = 'P2'
  return parameters


def test_ScalerBase():
  """Test the Base Scaler Class."""

  class SB_filler(ScalerBase):
    """Class to fill in abstract methods."""

    def update_for_minimisation(self, apm, curvatures=False):
      """Fill in abstract method."""

    def perform_error_optimisation(self):
      """Fill in abstract method."""

    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
      """Fill in abstract method."""

  scalerbase = SB_filler()
  assert scalerbase.experiments is None
  assert scalerbase.params is None
  assert scalerbase.reflection_table == []
  assert scalerbase.Ih_table is None
  assert scalerbase.initial_keys == []

  with pytest.raises(AssertionError):
    scalerbase.Ih_table = [1.0]
  # Set an Ih table with a mock that conforms to the id
  Ih_table = Mock()
  Ih_table.id_ = "IhTableBase"
  scalerbase.Ih_table = Ih_table
  assert scalerbase.Ih_table is Ih_table

  scalerbase_directory = dir(scalerbase)
  assert '_scaling_subset' in scalerbase_directory

  # Test setting space group.
  new_sg = "P 1"
  scalerbase.space_group = new_sg
  assert scalerbase.space_group == space_group(new_sg)
  new_sg = space_group("P 3")
  scalerbase.space_group = new_sg
  assert scalerbase.space_group is new_sg
  with pytest.raises(AssertionError):
    scalerbase.space_group = 1

  rt = flex.reflection_table()
  rt['inverse_scale_factor'] = flex.double([1.0])
  rt['inverse_scale_factor_variance'] = flex.double([1.0])
  rt['Ih_values'] = flex.double([1.0])
  rt['extra junk'] = flex.double([4.0])
  scalerbase._reflection_table = rt
  scalerbase.clean_reflection_table()
  assert not 'extra junk' in scalerbase.reflection_table
  assert 'inverse_scale_factor' in scalerbase.reflection_table
  assert 'inverse_scale_factor_variance' in scalerbase.reflection_table
  assert 'Ih_values' in scalerbase.reflection_table


def test_SingleScaler(test_reflections, test_experiments, test_params,
    mock_errormodel):
  """Test the single scaler class."""
  exp = create_scaling_model(test_params, test_experiments, test_reflections)
  singlescaler = SingleScalerBase(test_params, exp[0], test_reflections[0],
    scaled_id=2)

  # First test for things that are required upon initialisation.
  # Test that the var_cov matrix has been initialised to zero with the correct size
  assert singlescaler.var_cov_matrix.n_rows == 2
  assert singlescaler.var_cov_matrix.n_cols == 2
  assert singlescaler.var_cov_matrix.non_zeroes == 0
  assert singlescaler.space_group.info().symbol_and_number() == 'P 1 2 1 (No. 3)'
  assert singlescaler.consecutive_refinement_order == (
    singlescaler.experiments.scaling_model.consecutive_refinement_order)

  rt = singlescaler.reflection_table
  # Test that an id has been given.
  assert list(rt['id']) == [2] * rt.size()

  # Test that bad reflections are removed.
  assert list(rt.get_flags(rt.flags.integrated)) == [True, True, False, False]
  assert list(rt.get_flags(rt.flags.excluded_for_scaling)) == [
    False, False, True, True]
  assert list(rt.get_flags(rt.flags.outlier_in_scaling)) == [
    False, False, False, False]
  assert list(rt['inverse_scale_factor']) == [1.0] * rt.size()

  # Test for correct choice of intensities.
  new_rt = SingleScalerBase._select_optimal_intensities(rt, 'prf')
  assert list(new_rt['intensity']) == list(rt['intensity.prf.value'])
  assert list(new_rt['variance']) == list(rt['intensity.prf.variance'])
  new_rt = SingleScalerBase._select_optimal_intensities(rt, 'sum')
  assert list(new_rt['intensity']) == list(rt['intensity.sum.value'])
  assert list(new_rt['variance']) == list(rt['intensity.sum.variance'])
  # If bad choice, currently return the prf values.
  new_rt = SingleScalerBase._select_optimal_intensities(rt, 'bad')
  assert list(new_rt['intensity']) == list(rt['intensity.prf.value'])
  assert list(new_rt['variance']) == list(rt['intensity.prf.variance'])

  # Test that normalised Es are set (defer test of calculation to separate test)
  assert 'Esq' in rt
  assert 'intensity' in rt
  assert 'variance' in rt

  assert singlescaler.components == exp[0].scaling_model.components
  assert singlescaler.experiments == exp[0]
  assert singlescaler.params == test_params

  # Test configure_reflection_table?
  assert '_configure_reflection_table' in dir(singlescaler)
  # # Need a non-KB model as it doesnt do anything here

  # Now test public methods.
  # Test select reflections for scaling - should have selected two.
  assert singlescaler.Ih_table.size == 2
  assert list(singlescaler.Ih_table.asu_miller_index) == (
    list(flex.miller_index([(1, 0, 0), (0, 0, 1)])))
  assert singlescaler.components['scale'].n_refl == 2
  assert singlescaler.components['decay'].n_refl == 2

  # Test apply selection - here select just first.
  singlescaler.apply_selection_to_SFs(flex.bool([True, False, False, False]))
  assert singlescaler.components['scale'].n_refl == 1
  assert singlescaler.components['decay'].n_refl == 1
  singlescaler.apply_selection_to_SFs(flex.bool([True, True, False, False]))
  singlescaler.components['scale'].parameters = flex.double([1.1])

  rt = singlescaler.reflection_table
  assert list(rt['inverse_scale_factor']) == [1.0] * rt.size()
  assert 'inverse_scale_factor_variance' not in rt
  singlescaler.expand_scales_to_all_reflections()
  rt = singlescaler.reflection_table
  assert list(rt['inverse_scale_factor']) == [1.1, 1.1, 1.0, 1.0]
  # Variance set to zero if not calculated
  assert list(rt['inverse_scale_factor_variance']) == [0.0] * rt.size()

  # Test update var cov.
  apm = scaling_active_parameter_manager(singlescaler.components,
    ['scale', 'decay'])
  apm.var_cov_matrix = flex.double([1.0, 0.1, 0.1, 0.5])
  apm.var_cov_matrix.reshape(flex.grid(2, 2))
  singlescaler.update_var_cov(apm)
  # Should calculate the var cov of the two valid reflections.
  singlescaler.expand_scales_to_all_reflections(calc_cov=True)
  rt = singlescaler.reflection_table
  assert list(rt['inverse_scale_factor_variance']) == list(calc_sf_variances(
    singlescaler.components, singlescaler.var_cov_matrix)) + [0.0, 0.0]

  # Test update error model
  mock_em = mock_errormodel
  singlescaler.update_error_model(mock_em)
  assert mock_em.refined_parameters == (
    singlescaler.experiments.scaling_model.configdict['error_model_parameters'])

  singlescaler.apply_error_model_to_variances()
  assert singlescaler.reflection_table['variance'] == mock_em.update_variances()

  new_sg = "P 1"
  singlescaler.space_group = new_sg
  assert singlescaler.space_group == space_group(new_sg)

  with mock.patch('dials.algorithms.scaling.scaler.ScalingRestraints') as ScalingRestraints:
    _ = singlescaler.calculate_restraints(apm)
    assert ScalingRestraints.call_count == 1
    _ = singlescaler.compute_restraints_residuals_jacobian(apm)
    assert ScalingRestraints.call_count == 2

  # Test restraints calculation calls? But should these be moved elsewhere?
  # Test outlier rejection call? Defer to test in scaling utilities?

@pytest.fixture
def mock_errormodel():
  """A mock error model."""
  em = MagicMock()
  em.refined_parameters = [1.0, 0.1]
  em.update_variances.return_value = flex.double([1.0, 1.1, 1.0, 1.0])
  return em


def test_singlescaler_updateforminimisation(test_reflections,
    mock_exp, test_params):
  """Test the update_for_minimisation method of the singlescaler."""
  #First test the standard case.
  single_scaler = SingleScalerBase(test_params, mock_exp, test_reflections[0])
  apm = test_apm()
  assert list(single_scaler.Ih_table.inverse_scale_factors) == [1.0, 1.0]
  assert list(single_scaler.Ih_table.Ih_values) == [1.0, 10.0]
  single_scaler.update_for_minimisation(apm)
  #Should set new scale factors, and calculate Ih and weights.
  bf = basis_function(apm).return_basis()
  assert list(single_scaler.Ih_table.inverse_scale_factors) == list(bf[0])
  assert list(single_scaler.Ih_table.Ih_values) != [1.0, 10.0]
  assert list(single_scaler.Ih_table.Ih_values) == list(
    single_scaler.Ih_table.intensities / bf[0])
  for i in range(apm.derivatives.n_rows):
    for j in range(apm.derivatives.n_cols):
      assert apm.derivatives[i, j] == bf[1][i, j]
  assert apm.derivatives.non_zeroes == bf[1].non_zeroes

  # Now test for case when using free_Ih_table.
  test_params.scaling_options.use_free_set = True
  test_params.scaling_options.free_set_percentage = 50.0
  # Expect the first reflection in the Ih table and the second in the free table.
  single_scaler = SingleScalerBase(test_params, mock_exp, test_reflections[0])
  apm = test_apm()
  assert list(single_scaler.Ih_table.inverse_scale_factors) == [1.0]
  assert list(single_scaler.Ih_table.Ih_values) == [1.0]
  single_scaler.update_for_minimisation(apm)
  bf = basis_function(apm).return_basis()
  #Should set new scale factors, and calculate Ih and weights.
  assert list(single_scaler.Ih_table.inverse_scale_factors) == [bf[0][0]]
  assert list(single_scaler.Ih_table.Ih_values) != [1.0]
  assert list(single_scaler.Ih_table.Ih_values) == [1.0 / bf[0][0]]
  assert list(single_scaler.Ih_table.free_Ih_table.inverse_scale_factors) == [
    1.0 * bf[0][1]]
  assert list(single_scaler.Ih_table.free_Ih_table.Ih_values) != [1.0]
  assert list(single_scaler.Ih_table.free_Ih_table.Ih_values) == [10.0 / bf[0][1]]
  for i in range(apm.derivatives.n_rows):
    for j in range(apm.derivatives.n_cols):
      assert apm.derivatives[i, j] == bf[1][i, j]
  assert apm.derivatives.non_zeroes == 1

@pytest.fixture
def mock_scaling_component():
  """Mock scaling component to allow creation of a scaling model."""
  component = MagicMock()
  component.n_params = 2
  component.inverse_scales = flex.double([0.9, 1.1])
  component.derivatives = sparse.matrix(2, 2)
  component.derivatives[0, 0] = 0.5
  component.derivatives[1, 0] = 0.4
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
  exp.scaling_model.configure_reflection_table.side_effect = side_effect_config_table
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  exp.crystal = Crystal.from_dict(exp_dict)
  return exp

@pytest.fixture
def mock_explist_2(mock_exp):
  """A mock experimentlist."""
  mock_explist = []
  mock_explist.append(mock_exp)
  mock_explist.append(mock_exp)
  return mock_explist


class test_apm(object):
  """Mock-up of an apm for testing."""
  def __init__(self):
    self.var_list = [2.0]
    self.var_cov_matrix = flex.double(self.var_list)
    self.var_cov_matrix.reshape(flex.grid(1, 1))
    self.n_active_params = 1
    self.curvatures = []
    self.derivatives = sparse.matrix(1, 1)
    self.n_obs = 2
    self.constant_g_values = flex.double(2, 1.0)
    self.components_list = ['scale']
    self.components = {'scale': {'object' : mock_scaling_component(),
      'n_params' : 1, 'start_idx' : 0}}

def test_scaler_update_var_cov(test_reflections, mock_exp, test_params):
  """Test the update variance covariance matrix."""
  single_scaler = SingleScalerBase(test_params, mock_exp, test_reflections[0])
  #single_scaler.var_cov_matrix = sparse.matrix(2, 2)
  apm = Mock()
  var_list = [1.0, 0.1, 0.1, 0.5]
  apm.var_cov_matrix = flex.double(var_list)
  apm.var_cov_matrix.reshape(flex.grid(2, 2))
  apm.n_active_params = 2
  single_scaler.update_var_cov(apm)
  assert single_scaler.var_cov_matrix[0, 0] == var_list[0]
  assert single_scaler.var_cov_matrix[0, 1] == var_list[1]
  assert single_scaler.var_cov_matrix[1, 0] == var_list[2]
  assert single_scaler.var_cov_matrix[1, 1] == var_list[3]
  assert single_scaler.var_cov_matrix.non_zeroes == 4

  #Repeat but only with setting a subset - should only update one element
  single_scaler = SingleScalerBase(test_params, mock_exp, test_reflections[0])
  apm = test_apm()
  single_scaler.update_var_cov(apm)
  assert single_scaler.var_cov_matrix.non_zeroes == 1
  assert single_scaler.var_cov_matrix[0, 0] == apm.var_list[0]
  assert single_scaler.var_cov_matrix.n_cols == 2
  assert single_scaler.var_cov_matrix.n_rows == 2
  assert single_scaler.var_cov_matrix.non_zeroes == 1


@pytest.fixture
def mock_multiapm():
  """Create a mock-up of the multi apm for testing the update_for_minimisation
  of the multiscaler."""

  class mock_apm_class(object):
    """A mock apm class to hold a list."""
    def __init__(self, apm_list):
      self.apm_list = apm_list
      self.apm_data = {0 : {'start_idx' : 0}, 1: {'start_idx' : 2}}
      self.n_active_params = 4

  apm = mock_apm_class([Mock(), Mock()])
  return apm

@pytest.fixture
def mock_singlescaler(test_reflections_Ihtable, test_sg):
  """Mock singlescaler to use for testing multiscalers."""
  single_scaler = Mock()
  single_scaler.space_group = test_sg
  single_scaler.initial_keys = ['intensity', 'variance']
  single_scaler.reflection_table = test_reflections_Ihtable[0]
  single_scaler.Ih_table = SingleIhTable(test_reflections_Ihtable[0], test_sg)
  return single_scaler


def test_MultiScalerBase(mock_singlescaler, mock_explist_2, test_params,
    test_reflections_Ihtable):
  """Unit tests for the MultiScalerBase class."""

  class MSB_filler(MultiScalerBase):
    """Class to fill in abstract methods."""

    def join_multiple_datasets(self):
      """Fill in abstract method."""

    def update_for_minimisation(self, apm, curvatures=False):
      """Fill in abstract method."""

  singlescalers = [mock_singlescaler, mock_singlescaler]
  multiscaler = MSB_filler(test_params, mock_explist_2, singlescalers)

  assert multiscaler.single_scalers is singlescalers
  assert multiscaler.space_group is singlescalers[0].space_group
  assert multiscaler.active_scalers is None
  assert multiscaler.params is test_params
  assert multiscaler.experiments is mock_explist_2[0]
  assert isinstance(multiscaler.Ih_table, JointIhTable)

  # Set the scalers to active to allow calling of functions below
  multiscaler.active_scalers = multiscaler.single_scalers

  apm = Mock()
  fp = 'dials.algorithms.scaling.scaler.'
  with mock.patch(fp+'MultiScalingRestraints.calculate_restraints',
    return_value='test_restr') as ScalingRestraints:
    r = multiscaler.calculate_restraints(apm)
    assert ScalingRestraints.call_count == 1
    assert r == 'test_restr'
  with mock.patch(fp+'MultiScalingRestraints.calculate_jacobian_restraints',
    return_value='test_jacob') as ScalingRestraints:
    r = multiscaler.compute_restraints_residuals_jacobian(apm)
    assert ScalingRestraints.call_count == 1
    assert r == 'test_jacob'

  # Test calls for updating of individual error models.
  multiscaler.apply_error_model_to_variances()
  assert mock_singlescaler.apply_error_model_to_variances.call_count == 2

  # Test calls for individual reflection selection.
  multiscaler.select_reflections_for_scaling()
  assert mock_singlescaler.select_reflections_for_scaling.call_count == 2

  # Test calls for individual scale application to all reflections.
  multiscaler.expand_scales_to_all_reflections()
  assert mock_singlescaler.expand_scales_to_all_reflections.call_count == 2

  def outlier_rej_side_effect(*args):
    """Side effect for overriding the call to reject_outliers."""
    return args[0]

  with mock.patch(fp+'reject_outliers',
    side_effect=outlier_rej_side_effect) as rejout:
    rt = multiscaler.round_of_outlier_rejection(test_reflections_Ihtable[0])
    assert rt is test_reflections_Ihtable[0]
    assert rejout.call_count == 1
    assert rejout.call_args_list == [call(test_reflections_Ihtable[0],
      multiscaler.space_group, multiscaler.params.scaling_options.outlier_rejection,
      multiscaler.params.scaling_options.outlier_zmax)]

    multiscaler.join_datasets_from_scalers(multiscaler.single_scalers)
    assert rejout.call_count == 2
    expected_rt = flex.reflection_table()
    expected_rt.extend(mock_singlescaler.reflection_table)
    expected_rt.extend(mock_singlescaler.reflection_table)
    assert list(multiscaler.reflection_table) == list(expected_rt)

def basisfn_side_effect(*args):
  """Side effect to override various method calls."""
  #Expect it to return scales, derivatives and curvatures
  scales = flex.double([1.2, 1.4])
  derivatives = sparse.matrix(2, 2)
  derivatives[0, 0] = 1.0
  derivatives[1, 1] = 0.3
  return [scales, derivatives, None]

def do_nothing_side_effect(*args):
  """Side effect to override various method calls."""
  pass

def basisfn_side_effect_free(*args):
  """Side effect to override various method calls."""
  #Expect it to return scales, derivatives and curvatures
  scales = flex.double([1.2, 1.4])
  derivatives = sparse.matrix(2, 2)
  derivatives[0, 0] = 0.2
  derivatives[1, 1] = 0.4
  return [scales, derivatives, 'mock_curvs']



def test_MultiScaler(mock_singlescaler, mock_explist_2, test_params,
  mock_multiapm):
  """Unit tests for the multiscaler class."""

  singlescalers = [mock_singlescaler, mock_singlescaler]
  multiscaler = MultiScaler(test_params, mock_explist_2, singlescalers)

  assert multiscaler.id_ == 'multi'
  assert multiscaler.single_scalers == singlescalers
  assert multiscaler.active_scalers == singlescalers
  standard_Ih_intensities = multiscaler.Ih_table.intensities

  # Test update error model call to Ih_table

  with mock.patch.object(multiscaler.Ih_table, 'update_error_model',
    autospec=True, side_effect=do_nothing_side_effect) as update_em:
    multiscaler.update_error_model(Mock())
    assert update_em.call_count == 1

  # Test call to join multiple datasets - should call the method of
  # multiscalerbase with the single scalers.
  fp = 'dials.algorithms.scaling.scaler.'
  with mock.patch(fp+'MultiScalerBase.join_datasets_from_scalers',
    side_effect=do_nothing_side_effect) as join_data:
    multiscaler.join_multiple_datasets()
    assert join_data.call_args_list == [call(multiscaler.single_scalers)]

  # Now test update_for_minimisation method

  basis_output = basisfn_side_effect()

  with mock.patch(fp+'basis_function.return_basis',
    side_effect=basisfn_side_effect),\
    mock.patch.object(multiscaler.Ih_table, 'calc_Ih', autospec=True,
    side_effect=do_nothing_side_effect) as calc_Ih:
      multiscaler.update_for_minimisation(mock_multiapm)
      assert multiscaler.Ih_table.size == 4
      assert mock_multiapm.derivatives.n_rows == 4
      assert mock_multiapm.derivatives.n_cols == 4
      for scaler in multiscaler.single_scalers:
        assert scaler.Ih_table.inverse_scale_factors == basis_output[0]
      offset = 2
      for i in range(basis_output[1].n_rows):
        for j in range(basis_output[1].n_cols):
          assert mock_multiapm.derivatives[i, j] == basis_output[1][i, j]
          assert mock_multiapm.derivatives[i+offset, j+offset] == basis_output[1][i, j]
      assert mock_multiapm.derivatives.non_zeroes == 2 * basis_output[1].non_zeroes
      assert calc_Ih.call_count == 1

  # Now repeat the tests for the using free set option.
  test_params.scaling_options.use_free_set = True
  test_params.scaling_options.free_set_percentage = 50.0
  multiscaler = MultiScaler(test_params, mock_explist_2, singlescalers)
  assert multiscaler.Ih_table.size == 2
  assert multiscaler.Ih_table.free_Ih_table.size == 2
  assert list(multiscaler.Ih_table.intensities) == [
    standard_Ih_intensities[0], standard_Ih_intensities[2]]
  assert list(multiscaler.Ih_table.free_Ih_table.intensities) == [
    standard_Ih_intensities[1], standard_Ih_intensities[3]]

  # Test the udpate_for_minimisation. Expect to compose a derivatives of
  # n_refl x n_param, where n_refl is the number in the working set.



  basis_output = basisfn_side_effect_free()

  with mock.patch(fp+'basis_function.return_basis',
    side_effect=basisfn_side_effect_free),\
    mock.patch.object(multiscaler.Ih_table, 'calc_Ih', autospec=True,
    side_effect=do_nothing_side_effect) as calc_Ih:
      multiscaler.update_for_minimisation(mock_multiapm)
      assert multiscaler.Ih_table.size == 2
      assert multiscaler.Ih_table.free_Ih_table.size == 2
      assert mock_multiapm.derivatives.n_rows == 2
      assert mock_multiapm.derivatives.n_cols == 4
      for scaler in multiscaler.single_scalers:
        assert scaler.Ih_table.inverse_scale_factors == basis_output[0]
      offset_x = 1
      offset_y = 2
      for i in range(1):
        for j in range(2):
          assert mock_multiapm.derivatives[i, j] == basis_output[1][i, j]
          assert mock_multiapm.derivatives[i+offset_x, j+offset_y] == basis_output[1][i, j]
      assert mock_multiapm.derivatives.non_zeroes == 2
      assert calc_Ih.call_count == 1


def test_TargetScaler(mock_singlescaler, mock_explist_2, test_params,
  mock_multiapm):
  """Test for successful initialisation of TargetedScaler."""

  #singlescalers = [mock_singlescaler, mock_singlescaler]
  targetscaler = TargetScaler(test_params, mock_explist_2, [mock_singlescaler],
    [mock_singlescaler])

  assert targetscaler.unscaled_scalers == [mock_singlescaler]
  assert targetscaler.unscaled_scalers is targetscaler.active_scalers

  # Test the update_for_minimisation method.
  basis_output = basisfn_side_effect()
  fp = 'dials.algorithms.scaling.scaler.'
  with mock.patch(fp+'basis_function.return_basis',
    side_effect=basisfn_side_effect),\
    mock.patch.object(targetscaler.Ih_table, 'calc_Ih', autospec=True,
    side_effect=do_nothing_side_effect) as calc_Ih:
      targetscaler.update_for_minimisation(mock_multiapm)
      assert targetscaler.Ih_table.size == 2
      assert mock_multiapm.apm_list[0].derivatives.n_rows == 2
      assert mock_multiapm.apm_list[0].derivatives.n_cols == 2
      for scaler in targetscaler.unscaled_scalers:
        assert scaler.Ih_table.inverse_scale_factors == basis_output[0]
      for i in range(basis_output[1].n_rows):
        for j in range(basis_output[1].n_cols):
          assert mock_multiapm.apm_list[0].derivatives[i, j] == basis_output[1][i, j]
      assert mock_multiapm.apm_list[0].derivatives.non_zeroes == basis_output[1].non_zeroes
      assert calc_Ih.call_count == 0

  # For the target scaler - the free set selection applies to the individual
  # unscaled scalers, the target scaler is not changed in size as this is not
  # refined against.
  test_params.scaling_options.use_free_set = True
  test_params.scaling_options.free_set_percentage = 50.0
  targetscaler = TargetScaler(test_params, mock_explist_2, [mock_singlescaler],
    [mock_singlescaler])
  assert targetscaler.Ih_table.size == 2

  # Test the udpate_for_minimisation. Expect to compose a derivatives of
  # n_refl x n_param, where n_refl is the number in the working set.
  basis_output = basisfn_side_effect_free()

  with mock.patch(fp+'basis_function.return_basis',
    side_effect=basisfn_side_effect_free),\
    mock.patch.object(targetscaler.Ih_table, 'calc_Ih', autospec=True,
    side_effect=do_nothing_side_effect) as calc_Ih:
      targetscaler.update_for_minimisation(mock_multiapm)
      for k, scaler in enumerate(targetscaler.unscaled_scalers):
        assert scaler.Ih_table.size == 1
        assert scaler.Ih_table.free_Ih_table.size == 1
        assert mock_multiapm.apm_list[k].derivatives.n_rows == 1
        assert mock_multiapm.apm_list[k].derivatives.n_cols == 2
        assert scaler.Ih_table.inverse_scale_factors == basis_output[0][0]
        assert scaler.Ih_table.free_Ih_table.inverse_scale_factors == basis_output[0][1]
        for i in range(1):
          for j in range(2):
            assert mock_multiapm.apm_list[k].derivatives[i, j] == basis_output[1][i, j]
        assert mock_multiapm.apm_list[k].derivatives.non_zeroes == 1
      assert calc_Ih.call_count == 0

  fp = 'dials.algorithms.scaling.scaler.'
  with mock.patch(fp+'MultiScalerBase.join_datasets_from_scalers',
    side_effect=do_nothing_side_effect) as join_data:
    targetscaler.join_multiple_datasets()
    expected = []
    expected.extend(targetscaler.single_scalers)
    expected.extend(targetscaler.unscaled_scalers)
    assert join_data.call_args_list == [call(expected)]

    targetscaler.params.scaling_options.target_model = True
    targetscaler.join_multiple_datasets()
    assert join_data.call_args_list[1] == call(targetscaler.unscaled_scalers)

def test_sf_variance_calculation(test_experiments, test_params):
  """Test the calculation of scale factor variances."""
  assert len(test_experiments) == 1
  experiments = create_scaling_model(test_params, test_experiments, [None])
  components = experiments[0].scaling_model.components
  rt = flex.reflection_table()
  d1 = 1.0
  d2 = 2.0
  d3 = 3.0
  rt['d'] = flex.double([d1, d2, d3])
  components['scale'].update_reflection_data(rt)
  components['scale'].calculate_scales_and_derivatives()
  assert list(components['scale'].derivatives.col(0)) == [
    (0, 1.0), (1, 1.0), (2, 1.0)]
  components['decay'].update_reflection_data(rt)
  components['decay'].calculate_scales_and_derivatives()
  assert list(components['decay'].derivatives.col(0)) == [(0, 1.0/(2.0*d1*d1)),
    (1, 1.0/(2.0*d2*d2)), (2, 1.0/(2.0*d3*d3))]
  var_cov = sparse.matrix(2, 2)
  a = 0.2
  b = 0.3
  c = 0.1
  var_cov[0, 0] = a
  var_cov[0, 1] = c
  var_cov[1, 0] = c
  var_cov[1, 1] = b
  variances = calc_sf_variances(components, var_cov)
  assert approx_equal(list(variances), [b/(4.0*(d1**4.0)) + c/((d1**2.0)) + a,
    b/(4.0*(d2**4.0)) + c/((d2**2.0)) + a,
    b/(4.0*(d3**4.0)) + c/((d3**2.0)) + a])
