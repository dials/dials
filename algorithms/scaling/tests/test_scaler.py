'''
This code tests the data managers and active parameter managers.
'''
import copy as copy
import pytest
from scitbx import sparse
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.model.scaling_model_factory import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager, create_apm
#from dials.algorithms.scaling.active_parameter_managers import \
#  multi_active_parameter_manager
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  calc_sf_variances, ScalerBase
from dials.algorithms.scaling.Ih_table import JointIhTable

def generated_refl():
  '''function to generate input for datamanagers'''
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


@pytest.fixture(scope='module')
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
  assert 'get_basis_function' in scalerbase_directory
  assert '_scaling_subset' in scalerbase_directory
  assert '_update_weights_for_scaling' in scalerbase_directory

  # Test map indices to asu.
  rt = flex.reflection_table()
  rt['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (2, 2, 2)])
  exp = generated_exp()
  param = generated_param()
  sorted_rt = scalerbase._map_indices_to_asu(rt, exp[0], param)
  assert list(sorted_rt['miller_index']) == list(flex.miller_index([(0, 0, 1),
    (1, 0, 0), (1, 0, 0), (2, 2, 2)]))

def test_SingleScaler():
  """Test the single scaler class."""
  refl, exp, param = generated_refl(), generated_exp(), generated_param()
  exp = create_scaling_model(param, exp, refl)
  singlescaler = SingleScalerBase(param, exp[0], refl[0], scaled_id=2)

  # First test for things that are required upon initialisation.
  # Test that the var_cov matrix has been initialised to zero with the correct size
  assert singlescaler.var_cov_matrix.n_rows == 2
  assert singlescaler.var_cov_matrix.n_cols == 2
  assert singlescaler.var_cov_matrix.non_zeroes == 0

  rt = singlescaler.reflection_table
  # Test that an id has been given.
  assert list(rt['id']) == [2] * rt.size()

  # Test that reflection table has been sorted by asu miller index.
  assert list(rt['asu_miller_index']) == list(flex.miller_index(
    [(0, 0, 1), (1, 0, 0), (2, 0, 0), (2, 2, 2)]))

  # Test that bad reflections are removed.
  assert list(rt.get_flags(rt.flags.integrated)) == [True, True, False, False]
  assert list(rt.get_flags(rt.flags.excluded_for_scaling)) == [
    False, False, True, True]
  assert list(rt.get_flags(rt.flags.outlier_in_scaling)) == [
    False, False, False, False]
  assert list(rt['inverse_scale_factor']) == [1.0] * rt.size()

  # Test for correct choice of intensities.
  new_rt = SingleScalerBase._select_optimal_intensities(rt, param)
  assert list(new_rt['intensity']) == list(rt['intensity.prf.value'])
  assert list(new_rt['variance']) == list(rt['intensity.prf.variance'])
  param.scaling_options.integration_method = 'sum'
  new_rt = SingleScalerBase._select_optimal_intensities(rt, param)
  assert list(new_rt['intensity']) == list(rt['intensity.sum.value'])
  assert list(new_rt['variance']) == list(rt['intensity.sum.variance'])

  # Test that normalised Es are set (defer test of calculation to separate test)
  assert 'Esq' in rt
  assert 'intensity' in rt
  assert 'variance' in rt

  assert singlescaler.components == exp[0].scaling_model.components
  assert singlescaler.experiments == exp[0]
  assert singlescaler.params == param

  # Test configure_reflection_table?
  assert '_configure_reflection_table' in dir(singlescaler)
  # # Need a non-KB model as it doesnt do anything here

  # Now test public methods.
  # Test select reflections for scaling - should have selected two.
  assert singlescaler.Ih_table.size == 2
  assert list(singlescaler.Ih_table.asu_miller_index) == (
    list(flex.miller_index([(0, 0, 1), (1, 0, 0)])))
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
  bf = singlescaler.get_basis_function(apm)
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

@pytest.fixture(scope='module')
def generated_KB_param():
  """Generate a phil scope for a kB model."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters

def generated_multi_input():
  """Generate a multiple dataset input for a multiscaler."""
  refl = generated_refl()
  exp = generated_exp(2)
  param = generated_param()
  refl.append(generated_refl()[0])
  return (refl, exp, param)


#@pytest.fixture(scope='module')
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

def test_MultiScaler():
  """Test the MultiScaler class."""

  refl, exp, param = generated_multi_input()
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
  old_Ih_values = list(scaler.Ih_table.Ih_values)
  scaler.update_for_minimisation(apm)
  # Check Ih tables were updated.
  new_I_0 = list(scaler.single_scalers[0].Ih_table.inverse_scale_factors)
  new_I_1 = list(scaler.single_scalers[1].Ih_table.inverse_scale_factors)
  scaler.single_scalers[0].update_for_minimisation(apm.apm_list[0])
  single_bf_0 = scaler.single_scalers[0].get_basis_function(apm.apm_list[0])
  assert new_I_0 == list(single_bf_0[0])
  scaler.single_scalers[1].update_for_minimisation(apm.apm_list[1])
  single_bf_1 = scaler.single_scalers[1].get_basis_function(apm.apm_list[1])
  assert new_I_1 == list(single_bf_1[0])
  # Check derivatives were correctly composed.
  for i in range(2):
    for j in range(2):
      assert apm.derivatives[i, j] == single_bf_0[1][i, j]
  for i in range(2):
    for j in range(2):
      assert apm.derivatives[i+2, j+2] == single_bf_1[1][i, j]
  assert apm.derivatives.non_zeroes == 8
  # Check Ih values were updated. # Defer testing of calculation to test_Ih
  assert list(scaler.Ih_table.Ih_values) != old_Ih_values

  # Test that expand_scales_to_all_reflections updates each dataset.
  apm.set_param_vals(flex.double([1.1, 0.0, 1.2, 0.0]))
  scaler.update_for_minimisation(apm)
  scaler.expand_scales_to_all_reflections()
  assert list(scaler.single_scalers[0].reflection_table[
    'inverse_scale_factor']) == [1.1, 1.1, 1.0, 1.0]
  assert list(scaler.single_scalers[1].reflection_table[
    'inverse_scale_factor']) == [1.2, 1.2, 1.0, 1.0]

  # Test join_multiple_datasets - joins then sorts by asu_miller_idx
  scaler.join_multiple_datasets()
  assert list(scaler.reflection_table['inverse_scale_factor']) == [
    1.1, 1.2, 1.1, 1.2, 1.0, 1.0, 1.0, 1.0]

  # Other methods to test - update_error_model, calc_merging_stats.


'''def calculate_jacobian_fd(target):
  """Calculate jacobian matrix with finite difference approach."""
  delta = 1.0e-6
  #apm = target.apm
  jacobian = sparse.matrix(target.get_num_matches(), target.apm.n_active_params)
  #iterate over parameters, varying one at a time and calculating the residuals
  for i in range(target.apm.n_active_params):
    target.apm.x[i] -= 0.5 * delta
    target.predict()
    R_low = (target.calculate_residuals()/target.weights)**0.5 #unweighted unsquared residual
    target.apm.x[i] += delta
    target.predict()
    R_upper = (target.calculate_residuals()/target.weights)**0.5 #unweighted unsquared residual
    target.apm.x[i] -= 0.5 * delta
    target.predict()
    fin_difference = (R_upper - R_low) / delta
    for j in range(fin_difference.size()):
      jacobian[j, i] = fin_difference[j]
  return jacobian

def test_target_jacobian_calc():
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_param())
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)

  apm = scaling_active_parameter_manager(scaler.components, ['decay', 'scale'])

  target = ScalingTarget(scaler, apm)

  fd_jacobian = calculate_jacobian_fd(target)
  print(fd_jacobian)
  r, jacobian, w = target.compute_residuals_and_gradients()
  print(jacobian)
  print(list(w))
  for i in range(0, 3):
    for j in range(0, 2):
      assert approx_equal(jacobian[i, j], fd_jacobian[i, j])'''

def test_sf_variance_calculation(generated_KB_param):
  """Test the calculation of scale factor variances."""
  test_experiments, params = generated_exp(), generated_KB_param
  assert len(test_experiments) == 1
  experiments = create_scaling_model(params, test_experiments, [None])
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
