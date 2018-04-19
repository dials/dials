'''
This code tests the data managers and active parameter managers.
'''
import pytest
from scitbx import sparse
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.scaling_library import create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager, create_apm
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  calc_sf_variances, ScalerBase
from dials.algorithms.scaling.Ih_table import JointIhTable

@pytest.fixture
def test_reflections():
  """Make a test reflection table."""
  return generated_refl()

@pytest.fixture
def test_experiments():
  """Make a test experiments list"""
  return generated_exp()

@pytest.fixture
def two_test_experiments():
  """Make a test experiments list"""
  return generated_exp(2)

@pytest.fixture
def two_test_reflections():
  """Make a test experiments list"""
  refl = generated_refl()
  refl.append(generated_refl()[0])
  return refl

@pytest.fixture(scope='module')
def test_params():
  """Make a test param phil scope."""
  return generated_param()

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
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, False, False]),
    reflections.flags.integrated)
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
  goniometer_2 = Goniometer((1.0, 1.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  if n > 1:
    for _ in range(0, n-1):
      experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer_2,
        detector=detector, crystal=crystal))
  return experiments

def generated_param():
  """Generate a param phil scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters


def test_ScalerBase():
  """Test the Base Scaler Class."""

  class SB_filler(ScalerBase):
    """Class to fill in abstract methods."""

    def update_for_minimisation(self, apm, curvatures=False):
      """Fill in abstract method."""
      pass

    def perform_error_optimisation(self):
      """Fill in abstract method."""
      pass

    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
      """Fill in abstract method."""
      pass

  scalerbase = SB_filler()
  assert scalerbase.experiments is None
  assert scalerbase.params is None
  assert scalerbase.reflection_table == []
  assert scalerbase.Ih_table is None
  assert scalerbase.initial_keys == []

  with pytest.raises(AssertionError):
    scalerbase.Ih_table = [1.0]

  scalerbase_directory = dir(scalerbase)
  assert '_scaling_subset' in scalerbase_directory
  assert '_update_weights_for_scaling' in scalerbase_directory

  # Test map indices to asu.
  rt = flex.reflection_table()
  rt['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (2, 2, 2)])

def test_SingleScaler(test_reflections, test_experiments, test_params):
  """Test the single scaler class."""
  exp = create_scaling_model(test_params, test_experiments, test_reflections)
  singlescaler = SingleScalerBase(test_params, exp[0], test_reflections[0],
    scaled_id=2)

  # First test for things that are required upon initialisation.
  # Test that the var_cov matrix has been initialised to zero with the correct size
  assert singlescaler.var_cov_matrix.n_rows == 2
  assert singlescaler.var_cov_matrix.n_cols == 2
  assert singlescaler.var_cov_matrix.non_zeroes == 0

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
  new_rt = SingleScalerBase._select_optimal_intensities(rt, test_params)
  assert list(new_rt['intensity']) == list(rt['intensity.prf.value'])
  assert list(new_rt['variance']) == list(rt['intensity.prf.variance'])
  test_params.scaling_options.integration_method = 'sum'
  new_rt = SingleScalerBase._select_optimal_intensities(rt, test_params)
  assert list(new_rt['intensity']) == list(rt['intensity.sum.value'])
  assert list(new_rt['variance']) == list(rt['intensity.sum.variance'])

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

  # Test update for minimisation - this should call the basis function, set
  # derivatives and calc Ih .
  apm.set_param_vals([1.1, 0.0])
  assert apm.derivatives is None
  singlescaler.update_for_minimisation(apm)
  bf = basis_function(apm).return_basis()
  assert list(singlescaler.Ih_table.inverse_scale_factors) == list(bf[0])
  for i in range(apm.derivatives.n_rows):
    for j in range(apm.derivatives.n_cols):
      assert apm.derivatives[i, j] == bf[1][i, j]
  # Check that Ih values were updated.
  new_Ih = list(singlescaler.Ih_table.Ih_values)
  singlescaler.Ih_table.calc_Ih()
  assert new_Ih == list(singlescaler.Ih_table.Ih_values)

  # Test update error model
  singlescaler.update_error_model([1.0, 0.0])
  assert [1.0, 0.00] == singlescaler.experiments.scaling_model.configdict[
    'error_model_parameters']

  # Test restraints calculation calls? But should these be moved elsewhere?
  # Test outlier rejection call? Defer to test in scaling utilities?


def test_MultiScaler(two_test_reflections, two_test_experiments, test_params):
  """Test the MultiScaler class."""

  refl, exp, param = two_test_reflections, two_test_experiments, test_params
  exp = create_scaling_model(param, exp, refl)
  scaler = create_scaler(param, exp, refl)
  assert scaler.id_ == 'multi'
  assert isinstance(scaler.Ih_table, JointIhTable)
  assert len(scaler.single_scalers) == 2
  for single_scaler in scaler.single_scalers:
    assert isinstance(single_scaler, SingleScalerBase)

  # Test updating for minimisation
  apm_factory = create_apm(scaler)
  apm = apm_factory.make_next_apm()
  apm.set_param_vals(flex.double([1.1, 0.1, 1.2, 0.2]))
  # Check individual Ih tables were updated and derivatives matrix correctly composed
  #old_Ih_values = list(scaler.Ih_table.Ih_values)
  scaler.update_for_minimisation(apm)
  # Check Ih tables were updated.
  new_I_0 = list(scaler.single_scalers[0].Ih_table.inverse_scale_factors)
  new_I_1 = list(scaler.single_scalers[1].Ih_table.inverse_scale_factors)
  single_bf_0 = basis_function(apm.apm_list[0]).return_basis()
  assert new_I_0 == list(single_bf_0[0])
  single_bf_1 = basis_function(apm.apm_list[1]).return_basis()
  assert new_I_1 == list(single_bf_1[0])
  # Check derivatives were correctly composed.
  for i in range(2):
    for j in range(2):
      assert apm.derivatives[i, j] == single_bf_0[1][i, j]
  for i in range(2):
    for j in range(2):
      assert apm.derivatives[i+2, j+2] == single_bf_1[1][i, j]
  assert apm.derivatives.non_zeroes == 8

  # Test that expand_scales_to_all_reflections updates each dataset.
  apm.set_param_vals(flex.double([1.1, 0.0, 1.2, 0.0]))
  scaler.update_for_minimisation(apm)
  scaler.expand_scales_to_all_reflections()
  assert list(scaler.single_scalers[0].reflection_table[
    'inverse_scale_factor']) == [1.1, 1.1, 1.0, 1.0]
  assert list(scaler.single_scalers[1].reflection_table[
    'inverse_scale_factor']) == [1.2, 1.2, 1.0, 1.0]

  # Test join_multiple_datasets
  scaler.join_multiple_datasets()
  assert list(scaler.reflection_table['inverse_scale_factor']) == [
    1.1, 1.1, 1.0, 1.0, 1.2, 1.2, 1.0, 1.0]

  # Other methods to test - update_error_model, calc_merging_stats.


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
