"""
Tests for scaling library module.
"""

import pytest
from libtbx import phil
from mock import Mock, patch
from dials.util.options import OptionParser
from dials.array_family import flex
from cctbx import miller, crystal
from cctbx.sgtbx import space_group
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.scaling_library import scale_single_dataset,\
  create_scaling_model, create_datastructures_for_structural_model
from dials.algorithms.scaling.model.model import KBScalingModel
from dials.algorithms.scaling.model.scaling_model_factory import \
  PhysicalSMFactory

@pytest.fixture
def test_reflections():
  """Make a test reflection table."""
  return generated_refl()

@pytest.fixture
def test_experiments():
  """Make a test experiments list"""
  return generated_exp()

@pytest.fixture()
def test_params():
  """Make a test param phil scope."""
  return generated_param()

@pytest.fixture()
def mock_cif():
  """Mock a cif file for testing loading data from a cif."""
  cif = Mock()
  cif.intensities = flex.double([1.0, 1.0])
  cif.indices = flex.miller_index([(1, 0, 0), (0, 0, 1)])
  cif.space_group = space_group("C 2y")
  return cif

def generated_refl():
  """Generate a reflection table."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 10.0, 1.0, 2.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 10.0, 1.0, 2.0])
  reflections['intensity.sum.value'] = flex.double([10.0, 100.0, 100.0, 10.0, 10.0])
  reflections['intensity.sum.variance'] = flex.double([10.0, 100.0, 100.0, 10.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (0, 0, 1), (0, 0, 2)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 0.8, 2.0, 1.2]) #don't change
  reflections['lp'] = flex.double(5, 1.0)
  reflections['partiality'] = flex.double(5, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 10.0), (0.0, 0.0, 7.5)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0), (0.0, 0.1, 1.0), (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool(5, True),
    reflections.flags.integrated)
  return reflections

def generated_exp(n=1):
  """Generate an experiment list with two experiments."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 10], oscillation=[0.0, 1.0])
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
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  parameters.parameterisation.absorption_term = False
  parameters.parameterisation.n_resolution_bins = 1 # to stop example dataset
  # being overparameterised for array model refinement.
  return parameters

@pytest.mark.parametrize('model', ['physical', 'array'])
def test_scale_single_dataset(test_reflections, test_experiments, test_params,
    model):
  """Test completion of scaling."""
  scaled_reflections = scale_single_dataset(test_reflections, test_experiments,
    test_params, model=model)
  assert 'inverse_scale_factor' in scaled_reflections
  assert 'inverse_scale_factor_variance' in scaled_reflections

def test_create_scaling_model():
  """Test the create scaling model function."""

  # Test that one can create the correct scaling model with the phil param.
  for m in ['physical', 'array', 'KB']:
    params = generated_param()
    exp = generated_exp()
    rt = generated_refl()
    params.__inject__('model', m)
    new_exp = create_scaling_model(params, exp, [rt])
    assert new_exp[0].scaling_model.id_ == m

  # If a scaling model already exists, then nothing else should happen.
  params = generated_param()
  exp = generated_exp()
  rt = generated_refl()
  exp[0].scaling_model = PhysicalSMFactory().create(params, exp[0], rt)
  old_scaling_model = exp[0].scaling_model
  params.__inject__('model', 'KB')
  new_exp = create_scaling_model(params, exp, [rt])
  new_scaling_model = new_exp[0].scaling_model
  assert new_scaling_model is old_scaling_model # Should not modify original.

  # Test multiple datasets, where one already has a scaling model.
  exp = generated_exp(3)
  params = generated_param()
  rt = generated_refl()
  rt_2 = generated_refl()
  rt_3 = generated_refl()
  exp[0].scaling_model = PhysicalSMFactory().create(params, exp[0], rt)
  params.__inject__('model', 'KB')
  new_exp = create_scaling_model(params, exp, [rt, rt_2, rt_3])
  assert new_exp[0].scaling_model is exp[0].scaling_model
  assert isinstance(new_exp[1].scaling_model, KBScalingModel)
  assert isinstance(new_exp[2].scaling_model, KBScalingModel)


def mock_intensity_array_from_cif_file(cif):
  """Mock cif-intensity converter to replace call in create_datastructures..."""
  miller_set = miller.set(crystal_symmetry=crystal.symmetry(
    space_group=cif.space_group), indices=cif.indices, anomalous_flag=True)
  idata = miller.array(miller_set, data=cif.intensities)
  return idata

@patch('dials.algorithms.scaling.scaling_library.intensity_array_from_cif_file',
  side_effect=mock_intensity_array_from_cif_file)
def test_get_intensities_from_cif(_, test_reflections, test_experiments, mock_cif):
  """Test the conversion of a cif file to reflections and experiments suitable
  for scaling."""
  exp, refl = create_datastructures_for_structural_model([test_reflections],
    test_experiments, mock_cif)
  assert list(refl['intensity']) == [5.5, 5.5, 5.5, 5.5]
  assert list(refl['miller_index']) == [(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (0, 0, 1)]
  assert exp.scaling_model.is_scaled is True
