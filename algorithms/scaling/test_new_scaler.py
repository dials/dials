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

def side_effect_update_var(variances, intensities):
    """Side effect to mock configure reflection table
    call during initialisation."""
    return flex.double(range(1, len(variances)+1))

@pytest.fixture
def mock_errormodel():
  """A mock error model."""
  em = MagicMock()
  em.refined_parameters = [1.0, 0.1]
  em.update_variances.side_effect = side_effect_update_var
  return em

@pytest.fixture
def mock_errormodel2():
  """A mock error model."""
  em = MagicMock()
  em.refined_parameters = [1.0, 0.1]
  em.update_variances.side_effect = side_effect_update_var
  #return_value = flex.double(range(2, 9))
  return em

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

def generated_refl():
  """Create a reflection table suitable for splitting into blocks."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 2.0, 3.0, 4.0, 500.0, 6.0, 2.0, 2.0])
  reflections['variance'] = flex.double(8, 1.0)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (2, 0, 0), (0, 0, 1),
    (2, 2, 2), (1, 0, 0), (2, 0, 0), (1, 0, 0), (1, 0, 0)])
  reflections['d'] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6, 2.5, 2.5, 2.5])
  reflections['partiality'] = flex.double(8, 1.0)
  reflections['Esq'] = flex.double(8, 1.0)
  reflections['inverse_scale_factor'] = flex.double(8, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 8.0), (0.0, 0.0, 10.0), (0.0, 0.0, 12.0),
    (0.0, 0.0, 15.0), (0.0, 0.0, 15.0), (0.0, 0.0, 15.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0), (0.0, 0.1, 1.0),(0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool(8, True), reflections.flags.integrated)
  reflections.set_flags(flex.bool([False]*5 + [True] + [False]*2), reflections.flags.bad_for_scaling)
  reflections['id'] = flex.int(8, 0)
  reflections.experiment_identifiers()[0] = '0'
  return reflections

def generated_refl_2(exclude_refl=False):
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

def generated_refl_for_comb():
  """Create a reflection table suitable for splitting into blocks."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 2.0, 3.0, 4.0, 500.0, 6.0, 2.0, 2.0])
  reflections['variance'] = flex.double(8, 1.0)
  reflections['intensity.prf.value'] = flex.double([1.0, 3.0, 3.0, 4.0, 50.0, 6.0, 3.0, 2.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0])
  reflections['intensity.sum.value'] = flex.double([1.0, 4.0, 3.0, 4.0, 500.0, 6.0, 6.0, 2.0])
  reflections['intensity.sum.variance'] = flex.double(8, 1.0)
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (2, 0, 0), (0, 0, 1),
    (2, 2, 2), (1, 0, 0), (2, 0, 0), (1, 0, 0), (1, 0, 0)])
  reflections['d'] = flex.double([0.8, 2.1, 2.0, 1.4, 1.6, 2.5, 2.5, 2.5])
  reflections['partiality'] = flex.double(8, 1.0)
  reflections['Esq'] = flex.double(8, 1.0)
  reflections['inverse_scale_factor'] = flex.double(8, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 8.0), (0.0, 0.0, 10.0), (0.0, 0.0, 12.0),
    (0.0, 0.0, 15.0), (0.0, 0.0, 15.0), (0.0, 0.0, 15.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0), (0.0, 0.1, 1.0),(0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool(8, True), reflections.flags.integrated)
  reflections.set_flags(flex.bool([False]*5 + [True] + [False]*2), reflections.flags.bad_for_scaling)
  reflections['id'] = flex.int(8, 0)
  reflections.experiment_identifiers()[0] = '0'
  return reflections

def test_SingleScaler():
  p, e, r = (generated_param(), generated_exp(), generated_refl())
  exp = create_scaling_model(p, e, r)

  # test initialised correctly
  scaler = SingleScalerBase(p, exp[0], r)
  assert list(scaler.suitable_refl_for_scaling_sel) == [True] * 5 + [False] + [True] * 2
  # all 7 of the suitable should be within the scaling_subset
  assert list(scaler.scaling_subset_sel) == [True] * 7
  # one of these is not in the scaling selection due to being an outlier.
  assert list(scaler.scaling_selection) == [True] * 4 + [False] + [True] * 2
  assert list(scaler.outliers) == [False] * 4 + [True] + [False] * 2
  assert scaler.n_suitable_refl == 7

  # check for correct setup of global_Ih_table
  # block selection is order to extract out from suitable_reflections
  assert scaler.global_Ih_table.size == 7
  assert list(scaler.global_Ih_table.blocked_data_list[0].intensities) == \
    [3.0, 1.0, 500.0, 2.0, 2.0, 2.0, 4.0]
  block_selection = scaler.global_Ih_table.blocked_data_list[0].block_selections[0]
  assert list(block_selection) == [2, 0, 4, 5, 6, 1, 3]

  # check for correct setup of Ih_table
  assert scaler.Ih_table.size == 6
  assert list(scaler.Ih_table.blocked_data_list[0].intensities) == \
    [3.0, 1.0, 2.0, 2.0, 2.0, 4.0]
  block_selection = scaler.Ih_table.blocked_data_list[0].block_selections[0]
  assert list(block_selection) == [2, 0, 5, 6, 1, 3]

  # check for correct data/d_values in components
  d_suitable = r['d'].select(scaler.suitable_refl_for_scaling_sel)
  decay = scaler.experiments.scaling_model.components['decay']
  # first check 'data' contains all suitable reflections
  assert list(decay.data['d']) == list(d_suitable)
  # Now check 'd_values' (which will be used for minim.) matches Ih_table data
  assert list(decay.d_values[0]) == list(d_suitable.select(block_selection))

  # test expand to all reflections method. First check scales are all 1, then
  # update a component to simulate a minimisation result, then check that
  # scales are set only in all suitable reflections (as it may not be possible
  # to calculate scales for unsuitable reflections!)
  # Must also update the scales in the global_Ih_table
  assert list(scaler.reflection_table['inverse_scale_factor']) == [1.0] * 8
  scaler.experiments.scaling_model.components['scale'].parameters = flex.double([2.0])
  scaler.expand_scales_to_all_reflections(calc_cov=False)
  assert list(scaler.reflection_table['inverse_scale_factor']) == [2.0] * 5 + [1.0] + [2.0] * 2
  assert list(scaler.global_Ih_table.blocked_data_list[0].inverse_scale_factors) == [2.0] * 7

  # test make ready for scaling method
  # set some new outliers and check for updated datastructures
  scaler.outliers = flex.bool([False] * 3 + [True] * 2 + [False] * 2)
  scaler.make_ready_for_scaling(outlier=True)
  assert scaler.Ih_table.size == 5
  assert list(scaler.Ih_table.blocked_data_list[0].intensities) == \
    [3.0, 1.0, 2.0, 2.0, 2.0,]
  block_selection = scaler.Ih_table.blocked_data_list[0].block_selections[0]
  assert list(block_selection) == [2, 0, 5, 6, 1]
  assert list(decay.d_values[0]) == list(d_suitable.select(block_selection))

def generated_refl_2(exclude_refl=True):
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
  return reflections

def test_SingleScaler_update_for_minimisation():
  """Test the update_for_minimisation method of the singlescaler."""
  #test_params.scaling_options.nproc = 1
  p, e, r = (generated_param(), generated_exp(), generated_refl_2())
  exp = create_scaling_model(p, e, r)
  single_scaler = SingleScalerBase(p, exp[0], r)
  apm_fac = create_apm_factory(single_scaler)
  single_scaler.components['scale'].parameters /= 2.0
  apm = apm_fac.make_next_apm()

  Ih_table = single_scaler.Ih_table.blocked_data_list[0]
  Ih_table.calc_Ih()
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

@pytest.mark.xfail(reason='need to rework mcok error model')
def test_update_error_model(mock_errormodel, mock_errormodel2):
  """Test the update_error_model method"""
  p, e, r = (generated_param(), generated_exp(), generated_refl())
  exp = create_scaling_model(p, e, r)

  # test initialised correctly
  scaler = SingleScalerBase(p, exp[0], r)
  block = scaler.global_Ih_table.blocked_data_list[0]
  original_vars = block.variances
  # test update error model - should update weights in global Ih and vars in refl table
  scaler.update_error_model(mock_errormodel, apply_to_reflection_table=True)
  assert list(block.variances) == list(original_vars)
  newvars = flex.double(range(1, 8))
  assert list(block.block_selections[0]) == [2, 0, 4, 5, 6, 1, 3]
  # [1, 2, 3, 4, 5, 6, 7] < set these in ^ these positions (taking into account
  # the one non-suitable refl at index 5)
  new_vars_in_table = [2.0, 6.0, 1.0, 7.0, 3.0, 1.0, 4.0, 5.0]
  assert list(scaler.reflection_table['variance']) == list(new_vars_in_table)
  assert list(block.weights) == list(1.0/newvars)
  assert scaler.experiments.scaling_model.error_model is mock_errormodel

  # do again with second errormodel
  scaler.global_Ih_table.reset_error_model()
  scaler.update_error_model(mock_errormodel2, apply_to_reflection_table=True)
  assert list(block.variances) == list(original_vars)
  newvars = mock_errormodel2.update_variances()
  assert list(block.block_selections[0]) == [2, 0, 4, 5, 6, 1, 3]
  # [2, 3, 4, 5, 6, 7, 8] < set these in ^ these positions (taking into account
  # the one non-suitable refl at index 5)
  new_vars_in_table = [3.0, 7.0, 2.0, 8.0, 4.0, 1.0, 5.0, 6.0]
  assert list(scaler.reflection_table['variance']) == list(new_vars_in_table)
  assert list(block.weights) == list(1.0/newvars)
  assert scaler.experiments.scaling_model.error_model is mock_errormodel2

"""
  # test create Ih table method #implicity tested already

  # test update model data method #implicity tested already

  # test round of outlier rejection method #implicity tested already"""

def test_SingleScaler_combine_intensities():
  """test combine intensities method"""
  p, e, r = (generated_param(), generated_exp(), generated_refl_for_comb())
  exp = create_scaling_model(p, e, r)

  scaler = SingleScalerBase(p, exp[0], r)
  scaler.combine_intensities()

  # The input makes the profile intensities best - so check these are set in the
  # reflection table and global_Ih_table
  assert list(scaler.reflection_table['intensity']) == list(r['intensity.prf.value'])
  assert list(scaler.reflection_table['variance']) == list(r['intensity.prf.variance'])
  block = scaler.global_Ih_table.blocked_data_list[0]
  block_sel = block.block_selections[0]
  suitable = scaler.suitable_refl_for_scaling_sel
  assert list(block.intensities) == list(scaler.reflection_table['intensity'].select(
    suitable).select(block_sel))
  assert list(block.variances) == list(scaler.reflection_table['variance'].select(
    suitable).select(block_sel))
