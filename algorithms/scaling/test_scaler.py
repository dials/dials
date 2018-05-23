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
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager, create_apm_factory
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  calc_sf_variances, ScalerBase, MultiScalerBase, MultiScaler, TargetScaler
from dials.algorithms.scaling.Ih_table import IhTable

@pytest.fixture
def test_reflections():
  """Make a test reflection table."""
  return generated_refl()

@pytest.fixture
def test_2_reflections():
  """Make a test reflection table."""
  return [generated_refl_for_splitting_1()[0], generated_refl_for_splitting_2()[0]]

@pytest.fixture
def test_reflections_no_exclude():
  """Make a test reflection table."""
  return generated_refl_for_split()

@pytest.fixture
def test_reflections_Ihtable():
  """Make a test reflection table."""
  return generated_refl_Ihtable()

@pytest.fixture
def test_experiments():
  """Make a test experiments list"""
  return generated_exp()

@pytest.fixture
def test_2_experiments():
  """Make a test experiments list"""
  return generated_exp(2)

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
  reflections['d'] = flex.double([0.8, 2.0, 2.1, 0.1])
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  integrated_list = flex.bool([True, True, False, False])
  bad_list = flex.bool([False, False, True, True])
  reflections.set_flags(integrated_list, reflections.flags.integrated)
  reflections.set_flags(bad_list, reflections.flags.bad_for_scaling)
  return [reflections]

def generated_refl_for_splitting_1():
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
  reflections['intensity.prf.variance'] = flex.double(6, 1.0)
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
  reflections['intensity.prf.value'] = flex.double([7.0, 8.0, 9.0, 10.0, 11.0])
  reflections['intensity.prf.variance'] = flex.double(5, 1.0)
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

def generated_refl_for_split():
  """Generate a reflection table."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['intensity.sum.value'] = flex.double([10.0, 100.0, 1000.0, 10.0])
  reflections['intensity.sum.variance'] = flex.double([10.0, 100.0, 1000.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (2, 0, 0), (0, 0, 1),
    (2, 2, 2)]) #don't change
  #reflections['d'] = flex.double([0.8, 2.0, 2.0, 0.0]) #don't change
  reflections['d'] = flex.double([0.8, 2.1, 2.0, 0.1])
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 8.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  integrated_list = flex.bool(4, True)
  bad_list = flex.bool(4, False)
  reflections.set_flags(integrated_list, reflections.flags.integrated)
  reflections.set_flags(bad_list, reflections.flags.bad_for_scaling)
  return [reflections]


def generated_exp(n=1):
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
  if n > 1:
    for _ in range(n-1):
      experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
        detector=detector, crystal=crystal))
  return experiments

def generated_param():
  """Generate a param phil scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
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
  Ih_table.id_ = "IhTable"
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

def test_SingleScaler_splitintoblocks(test_reflections_no_exclude,
  test_experiments, test_params):
  test_params.model = 'physical'
  exp = create_scaling_model(test_params, test_experiments, test_reflections_no_exclude)
  test_params.scaling_options.n_proc = 2
  singlescaler = SingleScalerBase(test_params, exp[0], test_reflections_no_exclude[0],
    scaled_id=2)
  assert singlescaler.Ih_table.blocked_data_list
  assert list(singlescaler.components['decay'].d_values[0]) == [2.0, 0.8] #(#2 and #0)
  assert list(singlescaler.components['decay'].d_values[1]) == [2.1, 0.1] #(#1 and #3)
  shtt = singlescaler.components['absorption'].sph_harm_table.transpose()
  expected_harm1 = shtt.select_columns(flex.size_t([2, 0])).transpose()
  expected_harm2 = shtt.select_columns(flex.size_t([1, 3])).transpose()
  assert singlescaler.components['absorption'].harmonic_values[0] == expected_harm1
  assert singlescaler.components['absorption'].harmonic_values[1] == expected_harm2

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
  assert list(singlescaler.Ih_table.blocked_data_list[0].asu_miller_index) == (
    list(flex.miller_index([(0, 0, 1), (1, 0, 0)])))
  assert singlescaler.components['scale'].n_refl == [2]
  assert singlescaler.components['decay'].n_refl == [2]

  '''# Test apply selection - here select just first.
  singlescaler.apply_selection_to_SFs(flex.bool([True, False, False, False]))
  assert singlescaler.components['scale'].n_refl == [1]
  assert singlescaler.components['decay'].n_refl == [1]
  singlescaler.apply_selection_to_SFs(flex.bool([True, True, False, False]))'''
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
  Ih_table = single_scaler.Ih_table.blocked_data_list[0]
  assert list(Ih_table.inverse_scale_factors) == [1.0, 1.0]
  assert list(Ih_table.Ih_values) == [10.0, 1.0]
  single_scaler.update_for_minimisation(apm)
  #Should set new scale factors, and calculate Ih and weights.
  bf = basis_function(apm).calculate_scales_and_derivatives()
  assert list(Ih_table.inverse_scale_factors) == list(bf[0][0])
  assert list(Ih_table.Ih_values) != [1.0, 10.0]
  assert approx_equal(list(Ih_table.Ih_values), list(
    Ih_table.intensities / bf[0][0]))
  for i in range(Ih_table.derivatives.n_rows):
    for j in range(Ih_table.derivatives.n_cols):
      assert approx_equal(Ih_table.derivatives[i, j], bf[1][0][i, j])
  assert Ih_table.derivatives.non_zeroes == bf[1][0].non_zeroes

  # Don't test for free Ih_table until fixed.
  '''# Now test for case when using free_Ih_table.
  test_params.scaling_options.use_free_set = True
  test_params.scaling_options.free_set_percentage = 50.0
  # Expect the first reflection in the Ih table and the second in the free table.
  single_scaler = SingleScalerBase(test_params, mock_exp, test_reflections[0])
  apm = test_apm()
  assert list(single_scaler.Ih_table.inverse_scale_factors) == [1.0]
  assert list(single_scaler.Ih_table.Ih_values) == [1.0]
  single_scaler.update_for_minimisation(apm)
  bf = basis_function(apm).calculate_scales_and_derivatives()
  #Should set new scale factors, and calculate Ih and weights.
  assert list(single_scaler.Ih_table.inverse_scale_factors) == [bf[0][0]]
  assert list(single_scaler.Ih_table.Ih_values) != [1.0]
  assert list(single_scaler.Ih_table.Ih_values) == [1.0 / bf[0][0]]
  assert list(single_scaler.Ih_table.free_Ih_table.inverse_scale_factors) == [
    1.0 * bf[0][0][1]]
  assert list(single_scaler.Ih_table.free_Ih_table.Ih_values) != [1.0]
  assert list(single_scaler.Ih_table.free_Ih_table.Ih_values) == [10.0 / bf[0][0][1]]
  for i in range(single_scaler.Ih_table.derivatives.n_rows):
    for j in range(single_scaler.Ih_table.derivatives.n_cols):
      assert single_scaler.Ih_table.derivatives[i, j] == bf[1][0][i, j]
  assert single_scaler.Ih_table.derivatives.non_zeroes == 1'''

@pytest.fixture
def mock_scaling_component():
  """Mock scaling component to allow creation of a scaling model."""
  component = MagicMock()
  component.n_params = 2
  component.inverse_scales = [flex.double([0.9, 1.1])]
  component.derivatives = [sparse.matrix(2, 2)]
  component.derivatives[0][0, 0] = 0.5
  component.derivatives[0][1, 0] = 0.4
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
    self.derivatives = [sparse.matrix(1, 1)]
    self.n_obs = [2]
    self.constant_g_values = [flex.double(2, 1.0)]
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
def mock_singlescaler(test_reflections_Ihtable, test_sg):
  """Mock singlescaler to use for testing multiscalers."""
  single_scaler = Mock()
  single_scaler.space_group = test_sg
  single_scaler.initial_keys = ['intensity', 'variance']
  single_scaler.reflection_table = test_reflections_Ihtable[0]
  single_scaler.Ih_table = IhTable([(test_reflections_Ihtable[0], None)], test_sg)
  single_scaler.scaling_selection = flex.bool([True, True, False, False])
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

def do_nothing_side_effect(*args):
  """Side effect to override various method calls."""
  pass

def test_new_Multiscaler(test_2_reflections, test_2_experiments, test_params):
  """Test the setup of the Ih table and components for a multiscaler"""
  # Use the create_scaling_model and create_scaler helpers functions for ease.

  test_params.scaling_options.n_proc = 2
  test_params.model = 'physical'
  experiments = create_scaling_model(test_params, test_2_experiments,
    test_2_reflections)
  multiscaler = create_scaler(test_params, experiments, test_2_reflections)

  block_list = multiscaler.Ih_table.blocked_data_list
  block_sels = multiscaler.Ih_table.blocked_selection_list

  # Check that the expected miller indices have been sorted to the correct groups
  assert list(block_list[0].miller_index) == list(flex.miller_index([(0, 0, 1),
    (1, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 0)]))
  assert list(block_list[1].miller_index) == list(flex.miller_index([(2, 0, 0),
    (2, 0, 0), (2, 2, 2), (2, 0, 0), (2, 2, 2), (2, 2, 2)]))

  # For each block, matrices are ordered first by dataset id,
  # then by asu miller index
  expected_h_idx_1 = sparse.matrix(5, 2)
  expected_h_idx_1[0, 0] = 1
  expected_h_idx_1[1, 1] = 1
  expected_h_idx_1[2, 1] = 1
  expected_h_idx_1[3, 0] = 1
  expected_h_idx_1[4, 1] = 1

  expected_h_idx_2 = sparse.matrix(6, 2)
  expected_h_idx_2[0, 0] = 1
  expected_h_idx_2[1, 0] = 1
  expected_h_idx_2[2, 1] = 1
  expected_h_idx_2[3, 0] = 1
  expected_h_idx_2[4, 1] = 1
  expected_h_idx_2[5, 1] = 1

  assert block_list[0].h_index_matrix == expected_h_idx_1
  assert block_list[1].h_index_matrix == expected_h_idx_2

  #These are the selection lists to get the normalised values in the right order
  assert list(block_sels[0][0]) == [2, 0, 4] #dataset 0, block 0
  assert list(block_sels[0][1]) == [1, 5, 3] #dataset 0, block 1
  assert list(block_sels[1][0]) == [2, 4]    #dataset 1, block 0
  assert list(block_sels[1][1]) == [1, 0, 3] #dataset 1, block 1

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

  # Now test update_for_minimisation method. Make all
  # the parameters not one so that one can check the correct composition of
  # the inverse scales and derivatives by the method.
  apm_fac = create_apm_factory(multiscaler)
  multiscaler.single_scalers[0].components['scale'].parameters /= 2.0
  multiscaler.single_scalers[1].components['scale'].parameters *= 1.5
  apm = apm_fac.make_next_apm()
  multiscaler.update_for_minimisation(apm)
  # bf[0], bf[1] should be list of scales and derivatives
  bf1 = basis_function(apm.apm_list[0]).calculate_scales_and_derivatives()
  bf2 = basis_function(apm.apm_list[1]).calculate_scales_and_derivatives()
  expected_scales_for_block_1 = bf1[0][0]
  expected_scales_for_block_1.extend(bf2[0][0])
  expected_scales_for_block_2 = bf1[0][1]
  expected_scales_for_block_2.extend(bf2[0][1])

  expected_derivatives_for_block_1 = sparse.matrix(
    expected_scales_for_block_1.size(), apm.n_active_params)
  expected_derivatives_for_block_2 = sparse.matrix(
    expected_scales_for_block_2.size(), apm.n_active_params)

  expected_derivatives_for_block_1.assign_block(bf1[1][0], 0, 0)
  expected_derivatives_for_block_1.assign_block(bf2[1][0], bf1[1][0].n_rows,
    apm.apm_data[1]['start_idx'])
  expected_derivatives_for_block_2.assign_block(bf1[1][1], 0, 0)
  expected_derivatives_for_block_2.assign_block(bf2[1][1], bf1[1][1].n_rows,
    apm.apm_data[1]['start_idx'])


  block_list = multiscaler.Ih_table.blocked_data_list

  assert block_list[0].inverse_scale_factors == expected_scales_for_block_1
  assert block_list[1].inverse_scale_factors == expected_scales_for_block_2
  assert block_list[0].derivatives == expected_derivatives_for_block_1
  assert block_list[1].derivatives == expected_derivatives_for_block_2

def test_multiscaler_scaling(test_2_reflections, test_2_experiments, test_params):
  """Test the setup of the Ih table and components for a multiscaler.
  This should create some blocks with zero elements, but the algorithm should
  still complete."""
  # Use the create_scaling_model and create_scaler helpers functions for ease.
  test_2_reflections[1]['miller_index'][4] = flex.miller_index([(5,7,9)])[0]
  test_params.scaling_options.n_proc = 7
  test_params.scaling_refinery.engine = 'LevMar'
  # should split into 5 unique groups, but each dataset won't necessarily have
  # data in each block - the algorithm should still work!
  test_params.scaling_options.outlier_rejection = '0'
  test_params.model = 'KB'
  experiments = create_scaling_model(test_params, test_2_experiments,
    test_2_reflections)
  multiscaler = create_scaler(test_params, experiments, test_2_reflections)
  multiscaler.perform_scaling()

def test_new_TargetScaler(test_2_reflections, test_2_experiments, test_params):
  """Test the setup of the Ih table and components for a multiscaler"""
  # Use the create_scaling_model and create_scaler helpers functions for ease.

  test_params.scaling_options.n_proc = 2
  test_params.model = 'physical'
  experiments = create_scaling_model(test_params, test_2_experiments,
    test_2_reflections)
  experiments[0].scaling_model.set_scaling_model_as_scaled()
  target = create_scaler(test_params, experiments, test_2_reflections)
  assert isinstance(target, TargetScaler)


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
  assert list(components['scale'].derivatives[0].col(0)) == [
    (0, 1.0), (1, 1.0), (2, 1.0)]
  components['decay'].update_reflection_data(rt)
  components['decay'].calculate_scales_and_derivatives()
  assert list(components['decay'].derivatives[0].col(0)) == [(0, 1.0/(2.0*d1*d1)),
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
