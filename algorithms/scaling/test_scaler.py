'''
This code tests the data managers and active parameter managers.
'''
import pytest
import mock
from mock import Mock, MagicMock, call
from scitbx import sparse
from libtbx import phil
from libtbx.test_utils import approx_equal
from cctbx.sgtbx import space_group
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.array_family import flex
from dials.util.options import OptionParser
from dials.algorithms.scaling.scaling_library import create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager, create_apm_factory
from dials.algorithms.scaling.target_function import ScalingTargetFixedIH
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  calc_sf_variances, ScalerBase, MultiScalerBase, MultiScaler, TargetScaler,\
  NullScaler
from dials.algorithms.scaling.Ih_table import IhTable

@pytest.fixture
def test_reflections():
  """Make a test reflection table."""
  return generated_refl()

@pytest.fixture
def test_2_reflections():
  """Make a test reflection table."""
  return [generated_refl_for_splitting_1()[0],
    generated_refl_for_splitting_2()[0]]

@pytest.fixture
def test_reflections_no_exclude():
  """Make a test reflection table."""
  return generated_refl(exclude_refl=False)

@pytest.fixture
def test_reflections_Ihtable():
  """Make a test reflection table."""
  return generated_refl()

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

def refl_for_outlier_routine():
  """Generate a reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 1.0, 100.0, 1.0, 1.0])
  reflections['variance'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (1, 0, 0),
    (1, 0, 0), (1, 0, 0), (2, 0, 0)]) #don't change
  reflections['d'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0])
  reflections['Esq'] = flex.double([0.3, 1.0, 1.0, 1.0, 1.0])
  reflections['inverse_scale_factor'] = flex.double(5, 1.0)
  reflections['partiality'] = flex.double(5, 1.0)
  reflections['id'] = flex.int(5, 0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0)] * 5)
  integrated_list = flex.bool(5, True)
  bad_list = flex.bool(5, False)
  reflections.set_flags(integrated_list, reflections.flags.integrated)
  reflections.set_flags(bad_list, reflections.flags.bad_for_scaling)
  return reflections

def refl_for_error_optimisation():
  """Generate a reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double(range(0, 100))
  reflections['variance'] = flex.double(100, 1)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0)] * 100)
  reflections['d'] = flex.double(100, 2.0)
  Esq = flex.double(10, 0.1)
  Esq.extend(flex.double(90, 1.0))
  reflections['Esq'] = Esq
  reflections['inverse_scale_factor'] = flex.double(100, 1.0)
  reflections['id'] = flex.int(100, 0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, float(i))
    for i in range(0, 100)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0)] * 100)
  integrated_list = flex.bool(100, True)
  bad_list = flex.bool(100, False)
  reflections.set_flags(integrated_list, reflections.flags.integrated)
  reflections.set_flags(bad_list, reflections.flags.bad_for_scaling)
  return reflections

def generated_refl(exclude_refl=True):
  """Generate a reflection table."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['variance'] = flex.double([1.0, 10.0, 100.0, 1.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (2, 0, 0), (2, 2, 2)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0, 0.0]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.1, 0.1])
  reflections['Esq'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['inverse_scale_factor'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections['id'] = flex.int(4, 0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  if exclude_refl:
    integrated_list = flex.bool([True, True, False, False])
    bad_list = flex.bool([False, False, True, True])
  else:
    integrated_list = flex.bool(4, True)
    bad_list = flex.bool(4, False)
  reflections.set_flags(integrated_list, reflections.flags.integrated)
  reflections.set_flags(bad_list, reflections.flags.bad_for_scaling)
  return [reflections]

def generated_refl_for_splitting_1():
  """Create a reflection table suitable for splitting into blocks."""
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
  reflections['id'] = flex.int(6, 1)
  reflections.experiment_identifiers()[0] = '0'
  return [reflections]

def generated_refl_for_splitting_2():
  """Create another reflection table suitable for splitting into blocks."""
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
  reflections['id'] = flex.int(5, 2)
  reflections.experiment_identifiers()[1] = '1'
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
  experiments[0].identifier = '0'
  if n > 1:
    for i in range(n-1):
      experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
        detector=detector, crystal=crystal))
      experiments[i+1].identifier = str(i+1)
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

@pytest.fixture
def mock_apm():
  """mock parameter manager for testing var_cov_matrix setting."""
  apm = MagicMock()
  apm.var_cov_matrix = flex.double([2.0])
  apm.var_cov_matrix.reshape(flex.grid(1, 1))
  apm.n_active_params = 1
  apm.n_obs = [2]
  apm.curvatures = []
  apm.derivatives = [sparse.matrix(1, 1)]
  apm.components_list = ['scale']
  apm.components = {'scale': {'object' : mock_scaling_component(),
    'n_params' : 1, 'start_idx' : 0}}
  apm.constant_g_values = [flex.double(2, 1.0)]
  return apm

@pytest.fixture
def mock_errormodel():
  """A mock error model."""
  em = MagicMock()
  em.refined_parameters = [1.0, 0.1]
  em.update_variances.return_value = flex.double([1.1, 1.0, 0.1, 0.5])
  return em

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


def test_ScalerBase(test_params):
  """Test the Base Scaler Class."""

  class SB_filler(ScalerBase):
    """Class to fill in abstract methods."""

    def update_for_minimisation(self, apm, curvatures=False):
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

  # test scaling subset
  rt = generated_refl_for_splitting_1()[0]
  rt['intensity'] = rt['intensity.prf.value']
  rt['variance'] = rt['intensity.prf.variance']
  rt['Esq'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  test_params.reflection_selection.E2_range = 0.8, 5.0
  test_params.reflection_selection.d_range = 1.0, 5.0 #all but first
  test_params.reflection_selection.Isigma_range = 0.9, 5.5 #all but last
  sel = scalerbase._scaling_subset(rt, test_params)
  assert list(sel) == [False, True, True, True, True, False]
  rt['Esq'] = flex.double([1.0, 1.0, 1.0, 0.1, 6.0, 1.0])
  sel = scalerbase._scaling_subset(rt, test_params)
  assert list(sel) == [False, True, True, False, False, False]

  # test perform scaling in another test.
  rt1 = flex.reflection_table()
  scalerbase._reflection_table = rt1
  # test round_of_outlier_rejection - should just call reject outliers
  with mock.patch('dials.algorithms.scaling.scaler.reject_outliers') as \
    outlier_rej:
    scalerbase._params = test_params
    scalerbase.round_of_outlier_rejection()
    assert outlier_rej.call_count == 1
    outlier_rej.assert_called_with([rt1], scalerbase.space_group,
      scalerbase.params.scaling_options.outlier_rejection,
      scalerbase.params.scaling_options.outlier_zmax)


def test_SingleScaler(test_reflections, test_experiments, test_params,
    mock_errormodel, mock_apm):
  """Test the single scaler class."""
  test_params.scaling_options.nproc = 1
  exp = create_scaling_model(test_params, test_experiments, test_reflections)
  singlescaler = SingleScalerBase(test_params, exp[0], test_reflections[0])

  # First test for things that are required upon initialisation.
  # Test that the var_cov matrix has been initialised to zero with the correct size
  assert singlescaler.var_cov_matrix.n_rows == 2
  assert singlescaler.var_cov_matrix.n_cols == 2
  assert singlescaler.var_cov_matrix.non_zeroes == 0
  assert singlescaler.space_group.info().symbol_and_number() == 'P 1 2 1 (No. 3)'
  assert singlescaler.consecutive_refinement_order == (
    singlescaler.experiments.scaling_model.consecutive_refinement_order)

  rt = singlescaler.reflection_table

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
  # Test that Ih table was correctly set up by initialisation.
  assert singlescaler.Ih_table.size == 2
  assert list(singlescaler.Ih_table.blocked_data_list[0].asu_miller_index) == (
    list(flex.miller_index([(0, 0, 1), (1, 0, 0)])))
  assert singlescaler.components['scale'].n_refl == [2]
  assert singlescaler.components['decay'].n_refl == [2]

  singlescaler.components['scale'].parameters = flex.double([1.1])

  rt = singlescaler.reflection_table
  assert list(rt['inverse_scale_factor']) == [1.0] * rt.size()
  assert 'inverse_scale_factor_variance' not in rt

  # Test expand scales to all reflections
  singlescaler.expand_scales_to_all_reflections()
  rt = singlescaler.reflection_table
  assert list(rt['inverse_scale_factor']) == [1.1, 1.1, 1.0, 1.0]
  # Variance set to zero if not calculated
  assert list(rt['inverse_scale_factor_variance']) == [0.0] * rt.size()

  # Test update var cov.
  # First case - when var_cov_matrix is the full matrix.
  apm = Mock()
  apm.n_active_params = 2
  var_list = [1.0, 0.1, 0.1, 0.5]
  apm.var_cov_matrix = flex.double(var_list)
  apm.var_cov_matrix.reshape(flex.grid(2, 2))
  singlescaler.update_var_cov(apm)
  assert singlescaler.var_cov_matrix[0, 0] == var_list[0]
  assert singlescaler.var_cov_matrix[0, 1] == var_list[1]
  assert singlescaler.var_cov_matrix[1, 0] == var_list[2]
  assert singlescaler.var_cov_matrix[1, 1] == var_list[3]
  assert singlescaler.var_cov_matrix.non_zeroes == 4
  # Should calculate the var cov of the two valid reflections.
  singlescaler.expand_scales_to_all_reflections(calc_cov=True)
  rt = singlescaler.reflection_table
  assert list(rt['inverse_scale_factor_variance']) == pytest.approx(
    [1.5411376, 1.03695312, 0.0, 0.0]) #this was calculated previously by
    # a separate function, but behaviour changed so copied the verified result.

  # Second case - when var_cov_matrix is only part of full matrix.
  singlescaler = SingleScalerBase(test_params, exp[0], test_reflections[0])
  apm = scaling_active_parameter_manager(singlescaler.components, ['scale'])
  apm = mock_apm
  singlescaler.update_var_cov(apm)
  assert singlescaler.var_cov_matrix.non_zeroes == 1
  assert singlescaler.var_cov_matrix[0, 0] == 2.0
  assert singlescaler.var_cov_matrix.n_cols == 2
  assert singlescaler.var_cov_matrix.n_rows == 2
  assert singlescaler.var_cov_matrix.non_zeroes == 1
  singlescaler.expand_scales_to_all_reflections(calc_cov=True)
  assert list(rt['inverse_scale_factor_variance']) == [2.0, 2.0, 0.0, 0.0]

  #test create Ih table function with free set?

  # test reselect reflections for scaling - manually modify block selection list
  # and check components are updated
  singlescaler.Ih_table.blocked_selection_list[0] = [flex.size_t([0])]
  singlescaler.reselect_reflections_for_scaling()
  assert singlescaler.components['scale'].n_refl == [1]
  assert singlescaler.components['decay'].n_refl == [1]

  # Test update error model - should save to experiments object and call update
  # to Ih table
  with mock.patch.object(singlescaler.Ih_table, 'update_error_model',
    autospec=True) as update_em_1:
    mock_em = mock_errormodel
    singlescaler.update_error_model(mock_em)
    assert update_em_1.call_count == 1
  assert mock_em.refined_parameters == (
    singlescaler.experiments.scaling_model.configdict['error_model_parameters'])

  singlescaler.adjust_variances()
  assert singlescaler.reflection_table['variance'] == mock_em.update_variances()

  new_sg = "P 1"
  singlescaler.space_group = new_sg
  assert singlescaler.space_group == space_group(new_sg)

  rt1 = flex.reflection_table()
  rt1['inverse_scale_factor'] = flex.double([1.0])
  rt1['inverse_scale_factor_variance'] = flex.double([1.0])
  rt1['Ih_values'] = flex.double([1.0])
  rt1['Esq'] = flex.double([1.0])
  rt1['intensity'] = flex.double([1.0])
  rt1['variance'] = flex.double([1.0])
  rt1['extra junk'] = flex.double([4.0])
  singlescaler._reflection_table = rt1
  singlescaler.clean_reflection_tables()
  assert 'extra junk' not in singlescaler.reflection_table
  assert 'Esq' not in singlescaler.reflection_table
  assert 'inverse_scale_factor' in singlescaler.reflection_table
  assert 'inverse_scale_factor_variance' in singlescaler.reflection_table
  assert 'Ih_values' in singlescaler.reflection_table


def test_SingleScaler_error_optimisation(test_experiments, test_params):
  """Test perform error optimisation and error optimisation routine.
  The purpose is not to verify that the error model itself is correct, just
  that it is appropriately called and the scaler updated. Expect that the error
  model is minimised and updated in the Ih_table, experiments"""
  rt = refl_for_error_optimisation()
  test_params.scaling_options.nproc = 1
  exp = create_scaling_model(test_params, test_experiments, [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt)
  assert scaler.experiments.scaling_model.error_model is None
  initial_Ih_weights = list(scaler.Ih_table.blocked_data_list[0].weights)
  scaler.perform_error_optimisation()
  # Test that an error model has been set in the experiments
  assert scaler.experiments.scaling_model.error_model is not None
  # Test that the Ih_table weights have been updated.
  assert list(scaler.Ih_table.blocked_data_list[0].weights) != initial_Ih_weights

  # Now test the error_optimisation_routine
  rt = refl_for_error_optimisation()
  exp = create_scaling_model(test_params, test_experiments, [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt)

  with mock.patch.object(scaler, 'expand_scales_to_all_reflections'
    ) as expand_patch:
    with mock.patch.object(scaler, 'perform_error_optimisation'
      ) as optimise_patch:
      with mock.patch.object(scaler, 'reselect_reflections_for_scaling'
        ) as reselect_patch:
          scaler.error_optimisation_routine()
          assert expand_patch.call_count == 1
          assert optimise_patch.call_count == 1
          assert reselect_patch.call_count == 1
          scaler.error_optimisation_routine(make_ready_for_scaling=False)
          assert expand_patch.call_count == 2
          assert optimise_patch.call_count == 2
          assert reselect_patch.call_count == 1

def test_SingleScaler_outlier_rejection_routine(test_experiments, test_params):
  """Test outlier rejection routine. This function expects a scaler that
  has already been scaled - so first create this and check that outlier
  rejection is performed, a new Ih_table is created and then the right
  reflections reselected."""
  test_params.scaling_options.nproc = 1
  rt = refl_for_outlier_routine()
  test_params.scaling_options.outlier_rejection = 'standard'
  exp = create_scaling_model(test_params, test_experiments, [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt)
  # scaling_subset will select reflections [1, 2, 3, 4] then the rest will be
  # put into the Ih_table, and the indices of the blocked selection list
  # refers to selection from scaler.scaling_selection.
  assert list(scaler.Ih_table.blocked_selection_list[0][0]) == [0, 1, 2, 3]
  assert list(scaler.components['decay'].d_values[0]) == [2.0, 3.0, 4.0, 5.0]
  assert list(scaler.Ih_table.blocked_data_list[0].nonzero_weights) == [1, 2, 3, 4]
  # Now call outlier rejection routine - expect it to remove the reflection
  # in position 1.
  scaler.outlier_rejection_routine()
  assert list(scaler.Ih_table.blocked_data_list[0].nonzero_weights) == [1, 3, 4]
  assert list(scaler.Ih_table.blocked_selection_list[0][0]) == [0, 1, 2]
  assert list(scaler.components['decay'].d_values[0]) == [2.0, 4.0, 5.0]

  # Repeat but with make_ready_for_scaling=False
  rt = refl_for_outlier_routine() # create a fresh copy
  exp = create_scaling_model(test_params, test_experiments, [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt)
  assert list(scaler.Ih_table.blocked_data_list[0].nonzero_weights) == [1, 2, 3, 4]
  scaler.outlier_rejection_routine(make_ready_for_scaling=False)
  # Should still have found outlier
  assert list(scaler.reflection_table.get_flags(
    scaler.reflection_table.flags.outlier_in_scaling)) == [False, False,
      True, False, False]
  # Should not have created or updated Ih_table
  assert scaler.Ih_table == []
  #assert list(scaler.Ih_table.blocked_data_list[0].nonzero_weights) == [1, 2, 3, 4]

def test_SingleScaler_update_for_minimisation(test_reflections,
    test_experiments, test_params):
  """Test the update_for_minimisation method of the singlescaler."""
  test_params.scaling_options.nproc = 1
  exp = create_scaling_model(test_params, test_experiments, test_reflections[0])
  single_scaler = SingleScalerBase(test_params, exp[0], test_reflections[0])
  apm_fac = create_apm_factory(single_scaler)
  single_scaler.components['scale'].parameters /= 2.0
  apm = apm_fac.make_next_apm()

  Ih_table = single_scaler.Ih_table.blocked_data_list[0]
  assert list(Ih_table.inverse_scale_factors) == [1.0, 1.0]
  assert list(Ih_table.Ih_values) == [10.0, 1.0]
  single_scaler.update_for_minimisation(apm, 0)
  #Should set new scale factors, and calculate Ih and weights.
  bf = basis_function().calculate_scales_and_derivatives(apm.apm_list[0], 0)
  assert list(Ih_table.inverse_scale_factors) == list(bf[0])
  assert list(Ih_table.Ih_values) != [1.0, 10.0]
  assert approx_equal(list(Ih_table.Ih_values), list(
    Ih_table.intensities / bf[0]))
  for i in range(Ih_table.derivatives.n_rows):
    for j in range(Ih_table.derivatives.n_cols):
      assert approx_equal(Ih_table.derivatives[i, j], bf[1][i, j])
  assert Ih_table.derivatives.non_zeroes == bf[1].non_zeroes

def test_SingleScaler_split_into_blocks(test_reflections_no_exclude,
  test_experiments, test_params):
  """Test the scaler initialisation when nproc > 1 - the data in the Ih_Table
  and components should be correctly split into blocks."""
  test_params.model = 'physical'
  exp = create_scaling_model(test_params, test_experiments,
    test_reflections_no_exclude)
  test_params.scaling_options.nproc = 2
  singlescaler = SingleScalerBase(test_params, exp[0],
    test_reflections_no_exclude[0])
  assert singlescaler.Ih_table.blocked_data_list[0].size == 2
  assert singlescaler.Ih_table.blocked_data_list[1].size == 2
  assert list(singlescaler.components['decay'].d_values[0]) == [2.0, 0.8] #(#2 and #0)
  assert list(singlescaler.components['decay'].d_values[1]) == [2.1, 0.1] #(#1 and #3)
  shtt = singlescaler.components['absorption'].sph_harm_table#.transpose()
  expected_harm1 = shtt.select_columns(flex.size_t([2, 0])).transpose()
  expected_harm2 = shtt.select_columns(flex.size_t([1, 3])).transpose()
  assert singlescaler.components['absorption'].harmonic_values[0] == expected_harm1
  assert singlescaler.components['absorption'].harmonic_values[1] == expected_harm2

def test_MultiScalerBase(mock_singlescaler, mock_explist_2, test_params):
  """Unit tests for the MultiScalerBase class."""

  class MSB_filler(MultiScalerBase):
    """Class to fill in abstract methods."""

    def join_multiple_datasets(self):
      """Fill in abstract method."""

    def update_for_minimisation(self, apm, curvatures=False):
      """Fill in abstract method."""

  test_params.scaling_options.verbosity = 1

  singlescalers = [mock_singlescaler, mock_singlescaler]
  multiscaler = MSB_filler(test_params, mock_explist_2, singlescalers)

  assert multiscaler.single_scalers is singlescalers
  assert multiscaler.space_group is singlescalers[0].space_group
  assert multiscaler.active_scalers is None
  assert multiscaler.params is test_params
  assert multiscaler.experiments is mock_explist_2[0]

  # Set the scalers to active to allow calling of functions below
  multiscaler.active_scalers = multiscaler.single_scalers

  # Test calls for updating of individual error models.
  multiscaler.adjust_variances()
  assert mock_singlescaler.adjust_variances.call_count == 2

  # Test calls for individual reflection selection.
  #multiscaler.reselect_reflections_for_scaling()
  #assert mock_singlescaler.reselect_reflections_for_scaling.call_count == 2

  # Test calls for individual scale application to all reflections.
  multiscaler.expand_scales_to_all_reflections()
  assert mock_singlescaler.expand_scales_to_all_reflections.call_count == 2

  multiscaler.join_datasets_from_scalers(multiscaler.single_scalers)
  expected_rt = flex.reflection_table()
  expected_rt.extend(mock_singlescaler.reflection_table)
  expected_rt.extend(mock_singlescaler.reflection_table)
  assert list(multiscaler.reflection_table) == list(expected_rt)

def outlier_rej_side_effect(*args):
  """Side effect for overriding the call to reject_outliers."""
  return args[0]

def test_MultiScaler(test_2_reflections, test_2_experiments, test_params):
  """Test the setup of the Ih table and components for a multiscaler"""
  # Use the create_scaling_model and create_scaler helpers functions for ease.

  test_params.scaling_options.nproc = 2
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

  sscalers = multiscaler.active_scalers

  with mock.patch.object(multiscaler.Ih_table, 'update_error_model',
    autospec=True) as update_em_1:
      with mock.patch.object(sscalers[0].experiments.scaling_model, 'set_error_model',
        autospec=True) as update_em_2:
        with mock.patch.object(sscalers[1].experiments.scaling_model, 'set_error_model',
          autospec=True) as update_em_3:
          multiscaler.update_error_model(1.0)
          assert update_em_1.call_count == 1
          assert update_em_2.call_count == 1
          assert update_em_3.call_count == 1

  # Test call to join multiple datasets - should call the method of
  # multiscalerbase with the single scalers.
  fp = 'dials.algorithms.scaling.scaler.'
  with mock.patch(fp+'MultiScalerBase.join_datasets_from_scalers',
    ) as join_data:
    multiscaler.join_multiple_datasets()
    assert join_data.call_args_list == [call(multiscaler.single_scalers)]

  # Now test update_for_minimisation method. Make all
  # the parameters not one so that one can check the correct composition of
  # the inverse scales and derivatives by the method.
  apm_fac = create_apm_factory(multiscaler)
  multiscaler.single_scalers[0].components['scale'].parameters /= 2.0
  multiscaler.single_scalers[1].components['scale'].parameters *= 1.5
  apm = apm_fac.make_next_apm()
  multiscaler.update_for_minimisation(apm, 0)
  multiscaler.update_for_minimisation(apm, 1)
  # bf[0], bf[1] should be list of scales and derivatives
  s1, d1 = basis_function().calculate_scales_and_derivatives(apm.apm_list[0], 0)
  s2, d2 = basis_function().calculate_scales_and_derivatives(apm.apm_list[1], 0)
  s3, d3 = basis_function().calculate_scales_and_derivatives(apm.apm_list[0], 1)
  s4, d4 = basis_function().calculate_scales_and_derivatives(apm.apm_list[1], 1)
  expected_scales_for_block_1 = s1
  expected_scales_for_block_1.extend(s2)
  expected_scales_for_block_2 = s3
  expected_scales_for_block_2.extend(s4)

  expected_derivatives_for_block_1 = sparse.matrix(
    expected_scales_for_block_1.size(), apm.n_active_params)
  expected_derivatives_for_block_2 = sparse.matrix(
    expected_scales_for_block_2.size(), apm.n_active_params)

  expected_derivatives_for_block_1.assign_block(d1, 0, 0)
  expected_derivatives_for_block_1.assign_block(d2, d1.n_rows,
    apm.apm_data[1]['start_idx'])
  expected_derivatives_for_block_2.assign_block(d3, 0, 0)
  expected_derivatives_for_block_2.assign_block(d4, d3.n_rows,
    apm.apm_data[1]['start_idx'])

  block_list = multiscaler.Ih_table.blocked_data_list

  assert block_list[0].inverse_scale_factors == expected_scales_for_block_1
  assert block_list[1].inverse_scale_factors == expected_scales_for_block_2
  assert block_list[1].derivatives == expected_derivatives_for_block_2
  assert block_list[0].derivatives == expected_derivatives_for_block_1

  # Test the round_of_outlier_rejection method.
  assert multiscaler.reflection_table == []
  with mock.patch('dials.algorithms.scaling.scaler.reject_outliers',
    side_effect=outlier_rej_side_effect) as outlier_patch:
    multiscaler.round_of_outlier_rejection()
    assert outlier_patch.call_count == 1

def test_multiscaler_outlier_rejection_routine(test_params):
  """Test outlier rejection routine. This function expects a scaler that
  has already been scaled - so first create this and check that outlier
  rejection is performed, a new Ih_table is created and then the right
  reflections reselected."""
  test_params.scaling_options.nproc = 1
  rt = refl_for_outlier_routine()
  rt2 = refl_for_outlier_routine()
  rt2['intensity'][2] = 1.0
  test_params.scaling_options.outlier_rejection = 'standard'
  exp = create_scaling_model(test_params, generated_exp(), [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt, for_multi=True)
  exp = create_scaling_model(test_params, generated_exp(), [rt2])
  scaler2 = SingleScalerBase(test_params, exp[0], rt2, for_multi=True)
  multiscaler = MultiScaler(test_params, exp, [scaler, scaler2])
  # scaling_subset will select reflections [1, 2, 3, 4] then the rest will be
  # put into the Ih_table, and the indices of the blocked selection list
  # refers to selection from scaler.scaling_selection.
  assert list(multiscaler.Ih_table.blocked_selection_list[0][0]) == [0, 1, 2, 3]
  assert list(multiscaler.Ih_table.blocked_selection_list[1][0]) == [0, 1, 2, 3]
  assert list(multiscaler.single_scalers[0].components['decay'].d_values[0]
    ) == [2.0, 3.0, 4.0, 5.0]
  assert list(multiscaler.single_scalers[1].components['decay'].d_values[0]
    ) == [2.0, 3.0, 4.0, 5.0]
  assert list(multiscaler.Ih_table.blocked_data_list[0].nonzero_weights
    ) == [1, 2, 3, 4, 1, 2, 3, 4]
  assert multiscaler.Ih_table.blocked_data_list[0].size == 8
  # Now call outlier rejection routine - expect it to remove the reflection
  # in position 1.
  multiscaler.outlier_rejection_routine()
  assert list(multiscaler.Ih_table.blocked_data_list[0].nonzero_weights
    ) == [1, 3, 4, 1, 2, 3, 4]
  assert list(multiscaler.Ih_table.blocked_selection_list[0][0]) == [0, 1, 2]
  assert list(multiscaler.Ih_table.blocked_selection_list[1][0]) == [0, 1, 2, 3]
  assert list(multiscaler.single_scalers[0].components['decay'].d_values[0]
    ) == [2.0, 4.0, 5.0]
  assert list(multiscaler.single_scalers[1].components['decay'].d_values[0]
    ) == [2.0, 3.0, 4.0, 5.0]

  # Repeat but with make_ready_for_scaling=False - make a fresh copy of structures
  rt = refl_for_outlier_routine()
  rt2 = refl_for_outlier_routine()
  rt2['intensity'][2] = 1.0
  test_params.scaling_options.outlier_rejection = 'standard'
  test_params.scaling_options.nproc = 1
  exp = create_scaling_model(test_params, generated_exp(), [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt, for_multi=True)
  exp = create_scaling_model(test_params, generated_exp(), [rt2])
  scaler2 = SingleScalerBase(test_params, exp[0], rt2, for_multi=True)
  multiscaler = MultiScaler(test_params, exp, [scaler, scaler2])
  multiscaler.outlier_rejection_routine(make_ready_for_scaling=False)
  # Should still have found outlier
  assert list(multiscaler.single_scalers[0].reflection_table.get_flags(
    scaler.reflection_table.flags.outlier_in_scaling)) == [False, False,
      True, False, False]
  # Should not have created or updated Ih_table
  assert multiscaler.Ih_table == []
  #assert list(multiscaler.Ih_table.blocked_data_list[0].nonzero_weights
  ##  ) == [1, 2, 3, 1, 2, 3, 4, 4]

def test_multiscaler_scaling(test_2_reflections, test_2_experiments, test_params):
  """Test the setup of the Ih table and components for a multiscaler.
  This should create some blocks with zero elements, but the algorithm should
  still complete."""
  # Use the create_scaling_model and create_scaler helpers functions for ease.
  test_2_reflections[1]['miller_index'][4] = flex.miller_index([(5, 7, 9)])[0]
  test_params.scaling_options.nproc = 7
  test_params.scaling_refinery.engine = 'LevMar'
  # should split into 5 unique groups, but each dataset won't necessarily have
  # data in each block - the algorithm should still work!
  test_params.scaling_options.outlier_rejection = None
  test_params.model = 'KB'
  experiments = create_scaling_model(test_params, test_2_experiments,
    test_2_reflections)
  multiscaler = create_scaler(test_params, experiments, test_2_reflections)
  multiscaler.perform_scaling()

def test_new_TargetScaler(test_2_reflections, test_2_experiments, test_params):
  """Test the setup of the Ih table and components for a multiscaler"""
  # Use the create_scaling_model and create_scaler helpers functions for ease.

  test_params.scaling_options.nproc = 2
  test_params.model = 'physical'
  experiments = create_scaling_model(test_params, test_2_experiments,
    test_2_reflections)
  experiments[0].scaling_model.set_scaling_model_as_scaled()
  target = create_scaler(test_params, experiments, test_2_reflections)
  assert isinstance(target, TargetScaler)

  #Test join_multiple_datasets
  target.join_multiple_datasets()
  assert target.reflection_table.size() == \
    target.single_scalers[0].reflection_table.size() + \
    target.unscaled_scalers[0].reflection_table.size()

  #Test update for minimisation - should just call parent method with calc_Ih=False
  #Check that Ih's have not been updated
  mockapm = Mock()
  with mock.patch(
    'dials.algorithms.scaling.scaler.MultiScalerBase.update_for_minimisation'
    ) as update_patch:
    target.update_for_minimisation(mockapm)
    update_patch.assert_called_once_with(mockapm, False, calc_Ih=False)

  mockapm = Mock()
  with mock.patch('dials.algorithms.scaling.scaler.ScalerBase.perform_scaling'
    ) as scaling_patch:
    target.perform_scaling()
    scaling_patch.assert_called_once_with(target_type=ScalingTargetFixedIH,
      engine=None, max_iterations=None)

def test_targetscaler_outlier_rejection_routine(test_params):
  """Test outlier rejection routine. This function expects a scaler that
  has already been scaled - so first create this and check that outlier
  rejection is performed, a new Ih_table is created and then the right
  reflections reselected."""
  test_params.scaling_options.nproc = 1
  rt = refl_for_outlier_routine()
  rt2 = refl_for_outlier_routine()
  rt['intensity'][2] = 1.0
  test_params.scaling_options.outlier_rejection = 'standard'
  exp = create_scaling_model(test_params, generated_exp(), [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt, for_multi=True)
  exp = create_scaling_model(test_params, generated_exp(), [rt2])
  scaler2 = SingleScalerBase(test_params, exp[0], rt2, for_multi=True)
  targetscaler = TargetScaler(test_params, exp, [scaler], [scaler2])
  # scaling_subset will select reflections [1, 2, 3, 4] then the rest will be
  # put into the Ih_table, and the indices of the blocked selection list
  # refers to selection from scaler.scaling_selection.
  assert list(targetscaler.Ih_table.blocked_selection_list[0][0]) == [0, 1, 2, 3]
  assert list(targetscaler.unscaled_scalers[0].components['decay'].d_values[0]
    ) == [2.0, 3.0, 4.0, 5.0]
  assert list(targetscaler.Ih_table.blocked_data_list[0].nonzero_weights
    ) == [1, 2, 3, 4]
  # Now call outlier rejection routine - expect it to remove the reflection
  # in position 1.
  targetscaler.outlier_rejection_routine()
  assert list(targetscaler.Ih_table.blocked_data_list[0].nonzero_weights
    ) == [1, 3, 4]
  assert list(targetscaler.Ih_table.blocked_selection_list[0][0]) == [0, 1, 2]
  assert list(targetscaler.unscaled_scalers[0].components['decay'].d_values[0]
    ) == [2.0, 4.0, 5.0]

  # Repeat but with make_ready_for_scaling=False - make a fresh copy of structures
  rt = refl_for_outlier_routine()
  rt2 = refl_for_outlier_routine()
  rt['intensity'][2] = 1.0
  test_params.scaling_options.nproc = 1
  exp = create_scaling_model(test_params, generated_exp(), [rt])
  scaler = SingleScalerBase(test_params, exp[0], rt, for_multi=True)
  exp = create_scaling_model(test_params, generated_exp(), [rt2])
  scaler2 = SingleScalerBase(test_params, exp[0], rt2, for_multi=True)
  targetscaler = TargetScaler(test_params, exp, [scaler], [scaler2])
  targetscaler.outlier_rejection_routine(make_ready_for_scaling=False)
  # Should still have found outlier
  assert list(targetscaler.unscaled_scalers[0].reflection_table.get_flags(
    scaler.reflection_table.flags.outlier_in_scaling)) == [False, False,
      True, False, False]
  # Should not have created or updated Ih_table
  assert targetscaler.Ih_table == []
  #assert list(targetscaler.Ih_table.blocked_data_list[0].nonzero_weights
  #  ) == [1, 2, 3, 4]

def test_NullScaler(test_reflections, test_experiments, test_params):
  """Test for successful creation of NullScaler."""
  exp = create_scaling_model(test_params, test_experiments, test_reflections[0])
  _ = NullScaler(test_params, exp[0], test_reflections[0])
  # What exactly should be tested here?

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
  _, d = components['scale'].calculate_scales_and_derivatives()
  assert list(d.col(0)) == [
    (0, 1.0), (1, 1.0), (2, 1.0)]
  components['decay'].update_reflection_data(rt)
  s, d = components['decay'].calculate_scales_and_derivatives()
  assert list(d.col(0)) == [(0, 1.0/(2.0*d1*d1)),
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
