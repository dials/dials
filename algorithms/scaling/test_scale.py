"""
Test the command line script dials.scale, for successful completion.
Also tested are the command line files dials.plot_scaling_models and
dials.plot_scaling_outliers
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import cPickle as pickle
import pytest
from libtbx import easy_run, phil
from libtbx.utils import Sorry
from dxtbx.serialize import load
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.array_family import flex
from mock import Mock
import mock
from dials.util.options import OptionParser
from dials.command_line.scale import Script
from dials.algorithms.scaling.scaling_library import create_scaling_model

class run_one_scaling(object):
  """Class to run the dials.scale algorithm."""
  def __init__(self, pickle_path_list, sweep_path_list, extra_args, plot=True):
    args = ["dials.scale"] + pickle_path_list + sweep_path_list + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

    with open('scaled.pickle', 'rb') as fh:
      table = pickle.load(fh)

    assert 'inverse_scale_factor' in table
    assert 'inverse_scale_factor_variance' in table

    if plot:
      args = ["dials.plot_scaling_models", "scaled.pickle",
        "scaled_experiments.json"]
      command = " ".join(args)
      print(command)
      _ = easy_run.fully_buffered(command=command).raise_if_errors()
      args = ["dials.plot_scaling_outliers", "scaled.pickle",
        "scaled_experiments.json"]
      command = " ".join(args)
      print(command)
      _ = easy_run.fully_buffered(command=command).raise_if_errors()

def get_merging_stats(scaled_unmerged_mtz, anomalous=False, n_bins=20,
                  use_internal_variance=False, eliminate_sys_absent=False,
                  data_labels=None):
  import iotbx.merging_statistics
  i_obs = iotbx.merging_statistics.select_data(
    scaled_unmerged_mtz, data_labels=data_labels)
  i_obs = i_obs.customized_copy(anomalous_flag=False, info=i_obs.info())
  result = iotbx.merging_statistics.dataset_statistics(
    i_obs=i_obs,
    n_bins=n_bins,
    anomalous=anomalous,
    use_internal_variance=use_internal_variance,
    eliminate_sys_absent=eliminate_sys_absent,
  )
  return result

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
  """Generate the default scaling parameters object."""
  phil_scope = phil.parse('''
      include scope dials.command_line.scale.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  return parameters

def test_reflections():
  reflections = flex.reflection_table()
  reflections['intensity.sum.value'] = flex.double([1.0, 2.0, 3.0, 4.0])
  reflections['intensity.sum.variance'] = flex.double([1.0, 2.0, 3.0, 4.0])
  reflections['miller_index'] = flex.miller_index([(0, 0, 1), (0, 0, 1),
    (0, 0, 2), (0, 0, 2)])
  reflections['id'] = flex.int([0, 0, 0, 0])
  return reflections

def generate_test_input(n=1):
  reflections = []
  for _ in range(n):
    reflections.append(test_reflections())
  return generated_param(), generated_exp(n), reflections

def return_first_arg_side_effect(*args):
  """Side effect for overriding the call to reject_outliers."""
  return args[0]

def test_scale_merging_stats():
  """Test the merging stats method of dials.scale script"""
  params = generated_param()
  exp = generated_exp()
  reflections = flex.reflection_table()
  reflections['intensity.scale.value'] = flex.double([1.0, 2.0, 3.0, 4.0])
  reflections['intensity.scale.variance'] = flex.double([1.0, 2.0, 3.0, 4.0])
  reflections['miller_index'] = flex.miller_index([(0, 0, 1), (0, 0, 1),
    (0, 0, 2), (0, 0, 2)])
  reflections['id'] = flex.int([0, 0, 0, 0])
  reflections['inverse_scale_factor'] = flex.double([0.5, 0.4, 0.9, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0, 0, 0.5), (0, 0, 1.5),
    (0, 0, 2.5), (0, 0, 3.5)])
  reflections.set_flags(flex.bool(4, False), reflections.flags.bad_for_scaling)
  params.output.merging.nbins = 1
  script = Script(params, exp, [reflections])
  script.merging_stats()
  assert script.merging_statistics_result is not None

  # test for sensible return if small dataset with no equivalent reflections
  reflections['miller_index'] = flex.miller_index([(0, 0, 1), (0, 0, 2),
    (0, 0, 3), (0, 0, 4)])
  script = Script(params, exp, [reflections])
  script.merging_stats()
  assert script.merging_statistics_result is None

  #test expected behaviour of excluding methods
  #these methods require a scaling model, so make a simple Kb model
  reflections['miller_index'] = flex.miller_index([(0, 0, 1), (0, 0, 1),
    (0, 0, 2), (0, 0, 2)])
  with mock.patch('dials.command_line.scale.exclude_on_image_scale',
    side_effect=return_first_arg_side_effect) as \
    exclude_patch:
    params.output.exclude_on_image_scale = 0.6
    params.model = 'KB'
    script = Script(params, exp, [reflections])
    script.experiments = create_scaling_model(script.params, script.experiments,
      script.reflections)
    script.merging_stats()
    assert exclude_patch.call_count == 1
    assert exclude_patch.called_with(script.reflections, script.experiments,
      script.params.output.exclude_on_image_scale)
  with mock.patch('dials.command_line.scale.exclude_on_batch_rmerge',
    side_effect=return_first_arg_side_effect) as \
    exclude_patch:
    params.output.exclude_on_batch_rmerge = 2.0
    params.model = 'KB'
    script = Script(params, exp, [reflections])
    script.experiments = create_scaling_model(script.params, script.experiments,
      script.reflections)
    script.merging_stats()
    assert exclude_patch.call_count == 1
    assert exclude_patch.called_with(script.reflections, script.experiments,
      script.params.output.exclude_on_batch_rmerge)

def test_scale_script_prepare_input():
  """Test prepare_input method of scaling script."""

   #test the components of the scaling script directly with a test reflection
   #table, experiments list and params.

  params, exp, reflections = generate_test_input()
  #try to pass in unequal number of reflections and experiments
  reflections.append(test_reflections())
  script = Script(params, exp, reflections)
  with pytest.raises(Sorry):
    script.prepare_input()

  params, exp, reflections = generate_test_input()
  #Try to use use_datasets when not identifiers set
  params.dataset_selection.use_datasets = ['0']
  script = Script(params, exp, reflections)
  with pytest.raises(SystemExit):
    script.prepare_input()
  #Try to use use_datasets when not identifiers set
  params.dataset_selection.use_datasets = None
  params.dataset_selection.exclude_datasets = ['0']
  script = Script(params, exp, reflections)
  with pytest.raises(SystemExit):
    script.prepare_input()

  #Now make two experiments with identifiers and select on them
  params, exp, reflections = generate_test_input(n=2)
  exp[0].identifier = '0'
  reflections[0].experiment_identifiers()[0] = '0'
  exp[1].identifier = '1'
  reflections[1].experiment_identifiers()[0] = '1'
  assert reflections[0].are_experiment_identifiers_consistent([exp[0]])
  assert reflections[1].are_experiment_identifiers_consistent([exp[1]])
  params.dataset_selection.use_datasets = ['0']
  script = Script(params, exp, reflections)
  script.prepare_input()

  assert len(script.reflections) == 1
  assert script.reflections[0] is reflections[0]

  #Try again, this time excluding
  params, exp, reflections = generate_test_input(n=2)
  exp[0].identifier = '0'
  reflections[0].experiment_identifiers()[0] = '0'
  exp[1].identifier = '1'
  reflections[1].experiment_identifiers()[0] = '1'
  params.dataset_selection.exclude_datasets = ['0']
  script = Script(params, exp, reflections)
  script.prepare_input()

  assert len(script.reflections) == 1
  assert script.reflections[0] is reflections[1]

  #Try setting space group
  params, exp, reflections = generate_test_input(n=1)
  params.scaling_options.space_group = 'P1'
  script = Script(params, exp, reflections)
  script.prepare_input()
  assert script.experiments[0].crystal.get_space_group().type().number() == 1

  #Try having two unequal space groups
  params, exp, reflections = generate_test_input(n=2)
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " P 1"}
  crystal = Crystal.from_dict(exp_dict)
  exp[0].crystal = crystal
  script = Script(params, exp, reflections)
  with pytest.raises(Sorry):
    script.prepare_input()

  # Test cutting data
  params, exp, reflections = generate_test_input(n=1)
  reflections[0]['d'] = flex.double([2.5, 2.0, 1.0, 1.0])
  params.cut_data.d_min = 1.5
  script = Script(params, exp, reflections)
  script.prepare_input()
  r = script.reflections[0]
  assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
    False, False, True, True]
  params.cut_data.d_max = 2.25
  script = Script(params, exp, reflections)
  script.prepare_input()
  r = script.reflections[0]
  assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
    True, False, True, True]

  params, exp, reflections = generate_test_input(n=1)
  reflections[0]['partiality'] = flex.double([0.5, 0.8, 1.0, 1.0])
  params.cut_data.partiality_cutoff = 0.75
  script = Script(params, exp, reflections)
  script.prepare_input()
  r = script.reflections[0]
  assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [
    True, False, False, False]

@pytest.mark.skip(reason='Feature currently broken')
def test_scale_script_prepare_input_exclude_images():
  params, exp, reflections = generate_test_input(n=1)
  reflections[0]['xyzobs.px.value'] = flex.vec3_double([(0, 0, 0.5), (0, 0, 2.5),
    (0, 0, 4.5), (0, 0, 6.5)])
  params.cut_data.exclude_image_range = [0.0, 3.0]
  script = Script(params, exp, reflections)
  script.prepare_input()
  r = script.reflections[0]
  assert list(r.get_flags(r.flags.user_excluded_in_scaling)) == [True, True, False, False]

  params, exp, reflections = generate_test_input(n=2)
  params.cut_data.exclude_image_range = [0.0, 3.0]
  script = Script(params, exp, reflections)
  with pytest.raises(Sorry):
    script.prepare_input()

@pytest.mark.dataset_test
def test_scale_physical(dials_regression, tmpdir):
  """Test standard scaling of one dataset."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["model=physical", "merged_mtz=merged.mtz",
    "optimise_errors=False", "intensity_choice=profile", "unmerged_mtz=unmerged.mtz",
    "use_free_set=1", "outlier_rejection=simple"]

  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    assert os.path.exists("unmerged.mtz")
    assert os.path.exists("merged.mtz")

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats("unmerged.mtz")
    assert result.overall.r_pim < 0.024 #at 07/08/18, value was 0.0234517
    assert result.overall.cc_one_half > 0.9955 # at 07/08/18, value was 0.99597
    assert result.overall.n_obs > 2300 # at 07/08/18, was 2309

    # Try running again with the merged.mtz as a target, to trigger the
    # target_mtz option
    extra_args.append("target_mtz=merged.mtz")
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    result = get_merging_stats("unmerged.mtz")
    assert result.overall.r_pim < 0.024 #at 14/08/18, value was 0.023
    assert result.overall.cc_one_half > 0.9955 # at 14/08/18, value was 0.999
    assert result.overall.n_obs > 2100 # at 14/08/18, was 2123
    #FIXME in target_mtz, why are many more outliers rejected?

  # run again with the concurrent scaling option turned off and the 'standard'
  # outlier rejection
  extra_args = ["model=physical", "merged_mtz=merged.mtz",
    "unmerged_mtz=unmerged.mtz", "use_free_set=1", "outlier_rejection=standard",
    "concurrent=False", "intensity_choice=combine"]
  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    assert os.path.exists("scale_model.png")
    assert os.path.exists("absorption_surface.png")
    assert os.path.exists("outliers.png")

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats("unmerged.mtz")
    assert result.overall.r_pim < 0.035 #at 07/08/18, value was 0.034104
    assert result.overall.cc_one_half > 0.9935 # at 07/08/18, value was 0.99388
    assert result.overall.n_obs > 2310 # at 07/08/18, was 2319

@pytest.mark.dataset_test
def test_scale_optimise_errors(dials_regression, tmpdir):
  """Test standard scaling of one dataset with error optimisation."""
  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["model=physical", "optimise_errors=True"]
  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)

@pytest.mark.dataset_test
def test_scale_array(dials_regression, tmpdir):
  """Test a standard dataset - ideally needs a large dataset or full matrix
  round may fail. Currently turning off absorption term to avoid
  overparameterisation and failure of full matrix minimisation."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["model=array", "absorption_term=0", "full_matrix=0"]

  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    assert os.path.exists("decay_correction.png")
    assert os.path.exists("outliers.png")

@pytest.mark.dataset_test
def test_multi_scale(dials_regression, tmpdir):
  """Test standard scaling of two datasets."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path_1 = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path_1 = os.path.join(data_dir, "20_integrated_experiments.json")
  pickle_path_2 = os.path.join(data_dir, "25_integrated.pickle")
  sweep_path_2 = os.path.join(data_dir, "25_integrated_experiments.json")
  extra_args = ["unmerged_mtz=unmerged.mtz", "optimise_errors=False",
    "intensity_choice=profile", "outlier_rejection=simple"]

  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path_1, pickle_path_2],
      [sweep_path_1, sweep_path_2], extra_args)
    assert os.path.exists("scale_model_1.png")
    assert os.path.exists("scale_model_2.png")
    assert os.path.exists("absorption_surface_1.png")
    assert os.path.exists("absorption_surface_2.png")
    assert os.path.exists("outliers_1.png")
    assert os.path.exists("outliers_2.png")

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats("unmerged.mtz")
    expected_nobs = 5460
    assert abs(result.overall.n_obs - expected_nobs) < 10
    assert result.overall.r_pim < 0.022 #at 07/08/18, value was 0.021837
    assert result.overall.cc_one_half > 0.9975 # at 07/08/18, value was 0.99810

    #run again, optimising errors, and continuing from where last run left off.
    extra_args = ["optimise_errors=True", "unmerged_mtz=unmerged.mtz"]
    _ = run_one_scaling(["scaled.pickle"], ["scaled_experiments.json"],
      extra_args)
    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    # Note: error optimisation currently appears to give worse results here!
    result = get_merging_stats("unmerged.mtz")
    expected_nobs = 5520
    assert abs(result.overall.n_obs - expected_nobs) < 10
    assert result.overall.r_pim < 0.023 #at 07/08/18, value was 0.022722
    assert result.overall.cc_one_half > 0.9965 # at 07/08/18, value was 0.996925

@pytest.mark.dataset_test
def test_targeted_scaling(dials_regression, tmpdir):
  """Test the targeted scaling workflow."""
  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path_1 = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path_1 = os.path.join(data_dir, "20_integrated_experiments.json")
  pickle_path_2 = os.path.join(data_dir, "25_integrated.pickle")
  sweep_path_2 = os.path.join(data_dir, "25_integrated_experiments.json")

  extra_args = ["model=physical"]

  with tmpdir.as_cwd():

    args = ["dials.scale"] + [pickle_path_1] + [sweep_path_1] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

    experiments_list = load.experiment_list(
      "scaled_experiments.json", check_format=False)
    assert len(experiments_list.scaling_models()) == 1

    # Once individual has run, do targeted scaling of second dataset.
    # Use this as a chance to test the KB model as well.
    extra_args = ["model=KB"]
    args = ["dials.scale"] + ["scaled.pickle"] + ["scaled_experiments.json"] +\
      [pickle_path_2] + [sweep_path_2] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

    experiments_list = load.experiment_list(
      "scaled_experiments.json", check_format=False)
    assert len(experiments_list.scaling_models()) == 2
    assert experiments_list.scaling_models()[0].id_ == 'physical'
    assert experiments_list.scaling_models()[1].id_ == 'KB'

    extra_args = ["model=KB", "only_target=True"]
    args = ["dials.scale"] + ["scaled.pickle"] + ["scaled_experiments.json"] +\
      [pickle_path_2] + [sweep_path_2] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

@pytest.mark.dataset_test
def test_scale_cross_validate(dials_regression, tmpdir):
  """Test standard scaling of one dataset."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["cross_validation_mode=single", "nfolds=2", "full_matrix=0",
    "optimise_errors=0"]

  with tmpdir.as_cwd():
    args = ["dials.scale"] + [pickle_path] + [sweep_path] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()

  extra_args = ["cross_validation_mode=multi", "nfolds=2", "full_matrix=0",
    "optimise_errors=0", "parameter=absorption_term"]

  with tmpdir.as_cwd():
    args = ["dials.scale"] + [pickle_path] + [sweep_path] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()

  extra_args = ["cross_validation_mode=multi", "nfolds=2", "full_matrix=0",
    "optimise_errors=0", "parameter=model", 'parameter_values="physical array"']

  with tmpdir.as_cwd():
    args = ["dials.scale"] + [pickle_path] + [sweep_path] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()