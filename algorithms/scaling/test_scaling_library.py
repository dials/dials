"""
Tests for scaling library module.
"""
from __future__ import absolute_import, division, print_function

import pytest
from libtbx import phil
from mock import Mock, patch
from cctbx import miller, crystal
from cctbx.sgtbx import space_group
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.util.options import OptionParser
from dials.array_family import flex
from dials.algorithms.scaling.scaling_library import (
    scale_single_dataset,
    create_scaling_model,
    create_datastructures_for_structural_model,
    create_Ih_table,
    # calculate_merging_statistics,
    # calculate_single_merging_stats,
    choose_scaling_intensities,
    create_auto_scaling_model,
)
from dials.algorithms.scaling.model.model import KBScalingModel, PhysicalScalingModel
from dials.algorithms.scaling.model.scaling_model_factory import PhysicalSMFactory


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
    # these miller_idx/d_values don't make physical sense, but I didn't want to
    # have to write the tests for lots of reflections.
    reflections = flex.reflection_table()
    reflections["intensity.prf.value"] = flex.double([1.0, 10.0, 10.0, 1.0, 2.0])
    reflections["intensity.prf.variance"] = flex.double([1.0, 10.0, 10.0, 1.0, 2.0])
    reflections["intensity.sum.value"] = flex.double([10.0, 100.0, 100.0, 10.0, 10.0])
    reflections["intensity.sum.variance"] = flex.double(
        [10.0, 100.0, 100.0, 10.0, 10.0]
    )
    reflections["miller_index"] = flex.miller_index(
        [(1, 0, 0), (0, 0, 1), (1, 0, 0), (0, 0, 1), (0, 0, 2)]
    )  # don't change
    reflections["d"] = flex.double([0.8, 2.0, 0.8, 2.0, 1.2])  # don't change
    reflections["lp"] = flex.double(5, 1.0)
    reflections["partiality"] = flex.double([1.0, 1.0, 1.0, 1.0, 0.8])
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 5.0),
            (0.0, 0.0, 10.0),
            (0.0, 0.0, 10.0),
            (0.0, 0.0, 7.5),
        ]
    )
    reflections["s1"] = flex.vec3_double(
        [
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
            (0.0, 0.1, 1.0),
        ]
    )
    reflections.set_flags(flex.bool(5, True), reflections.flags.integrated)
    reflections["id"] = flex.int(5, 0)
    reflections.experiment_identifiers()[0] = str(0)
    return reflections


def generated_exp(n=1, scan=True, image_range=[0, 10]):
    """Generate an experiment list with two experiments."""
    experiments = ExperimentList()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 2.0],
        "space_group_hall_symbol": " C 2y",
    }
    crystal = Crystal.from_dict(exp_dict)
    if scan:
        scan = Scan(image_range=image_range, oscillation=[0.0, 1.0])
    else:
        scan = None
    beam = Beam(s0=(0.0, 0.0, 1.01))
    goniometer = Goniometer((1.0, 0.0, 0.0))
    goniometer_2 = Goniometer((1.0, 1.0, 0.0))
    detector = Detector()
    experiments.append(
        Experiment(
            beam=beam,
            scan=scan,
            goniometer=goniometer,
            detector=detector,
            crystal=crystal,
        )
    )
    experiments[0].identifier = "0"
    if n > 1:
        for i in range(0, n - 1):
            experiments.append(
                Experiment(
                    beam=beam,
                    scan=scan,
                    goniometer=goniometer_2,
                    detector=detector,
                    crystal=crystal,
                )
            )
            experiments[i + 1].identifier = str(i + 1)
    return experiments


def generated_param(absorption_term=False):
    """Generate a param phil scope."""
    phil_scope = phil.parse(
        """
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
  """,
        process_includes=True,
    )
    optionparser = OptionParser(phil=phil_scope, check_format=False)
    parameters, _ = optionparser.parse_args(
        args=[], quick_parse=True, show_diff_phil=False
    )
    parameters.parameterisation.absorption_term = absorption_term
    parameters.parameterisation.n_resolution_bins = 1  # to stop example dataset
    # being overparameterised for array model refinement.
    return parameters


@pytest.mark.parametrize("model", ["physical", "array"])
def test_scale_single_dataset(test_reflections, test_experiments, test_params, model):
    """Test completion of scaling."""
    scaled_reflections = scale_single_dataset(
        test_reflections, test_experiments, test_params, model=model
    )
    assert "inverse_scale_factor" in scaled_reflections
    assert "inverse_scale_factor_variance" in scaled_reflections
    # what about when no params supplied?


def test_scale_single_dataset_no_params_supplied(test_reflections, test_experiments):
    """Test when no params scope supplied."""
    scaled_reflections = scale_single_dataset(
        test_reflections, test_experiments, model="physical"
    )
    assert "inverse_scale_factor" in scaled_reflections
    assert "inverse_scale_factor_variance" in scaled_reflections


def test_create_scaling_model():
    """Test the create scaling model function."""

    # Test that one can create the correct scaling model with the phil param.
    for m in ["physical", "array", "KB"]:
        params = generated_param()
        exp = generated_exp()
        rt = generated_refl()
        params.__inject__("model", m)
        new_exp = create_scaling_model(params, exp, [rt])
        assert new_exp[0].scaling_model.id_ == m

    # If a scaling model already exists, then nothing else should happen.
    params = generated_param()
    exp = generated_exp()
    rt = generated_refl()
    exp[0].scaling_model = PhysicalSMFactory().create(params, exp[0], rt)
    old_scaling_model = exp[0].scaling_model
    params.__inject__("model", "KB")
    new_exp = create_scaling_model(params, exp, [rt])
    new_scaling_model = new_exp[0].scaling_model
    assert new_scaling_model is old_scaling_model  # Should not modify original.

    # Test multiple datasets, where one already has a scaling model.
    exp = generated_exp(3)
    params = generated_param()
    rt = generated_refl()
    rt_2 = generated_refl()
    rt_3 = generated_refl()
    exp[0].scaling_model = PhysicalSMFactory().create(params, exp[0], rt)
    params.__inject__("model", "KB")
    new_exp = create_scaling_model(params, exp, [rt, rt_2, rt_3])
    assert new_exp[0].scaling_model is exp[0].scaling_model
    assert isinstance(new_exp[1].scaling_model, KBScalingModel)
    assert isinstance(new_exp[2].scaling_model, KBScalingModel)

    # Now test overwrite_existing_models option
    params.overwrite_existing_models = True
    params.model = "physical"
    newer_exp = create_scaling_model(params, new_exp, [rt, rt_2, rt_3])
    for exp in newer_exp:
        assert isinstance(exp.scaling_model, PhysicalScalingModel)


def mock_intensity_array_from_cif_file(cif):
    """Mock cif-intensity converter to replace call in create_datastructures..."""
    miller_set = miller.set(
        crystal_symmetry=crystal.symmetry(space_group=cif.space_group),
        indices=cif.indices,
        anomalous_flag=True,
    )
    idata = miller.array(miller_set, data=cif.intensities)
    return idata


@patch(
    "dials.algorithms.scaling.scaling_library.intensity_array_from_cif_file",
    side_effect=mock_intensity_array_from_cif_file,
)
def test_get_intensities_from_cif(_, test_reflections, test_experiments, mock_cif):
    """Test the conversion of a cif file to reflections and experiments suitable
    for scaling."""
    exp, refl = create_datastructures_for_structural_model(
        [test_reflections], test_experiments, mock_cif
    )
    assert list(refl["intensity"]) == [5.5, 5.5, 5.5, 5.5]
    assert list(refl["miller_index"]) == [(1, 0, 0), (0, 0, 1), (1, 0, 0), (0, 0, 1)]
    assert exp.scaling_model.is_scaled is True


def return_len_refl_side_effect(*args):
    """Side effect for overriding the call to reject_outliers."""
    return args[0].size()


def return_data_side_effect(*args, **kwargs):
    """Side effect to return data from a miller array."""
    return kwargs["i_obs"].data()


def test_create_Ih_table(test_experiments, test_reflections):
    """Test the create_Ih_table function."""
    test_reflections["intensity"] = test_reflections["intensity.prf.value"]
    test_reflections["variance"] = test_reflections["intensity.prf.variance"]

    Ih_table = create_Ih_table(test_experiments, [test_reflections])

    # Test data has been sorted into one block as expected.
    assert list(Ih_table.blocked_data_list[0].asu_miller_index) == (
        [(0, 0, 1), (0, 0, 1), (0, 0, 2), (1, 0, 0), (1, 0, 0)]
    )

    selection = flex.bool([True, True, True, True, False])
    Ih_table = create_Ih_table(test_experiments, [test_reflections], [selection])
    assert list(Ih_table.blocked_data_list[0].asu_miller_index) == (
        [(0, 0, 1), (0, 0, 1), (1, 0, 0), (1, 0, 0)]
    )

    with pytest.raises(AssertionError):
        # Should fail if length of exp, refl and selection are different
        Ih_table = create_Ih_table(
            test_experiments, [test_reflections, test_reflections], [selection]
        )


def test_choose_scaling_intensities(test_reflections):
    """Test for correct choice of intensities."""
    test_refl = test_reflections
    intstr = "prf"
    new_rt = choose_scaling_intensities(test_refl, intstr)
    assert list(new_rt["intensity"]) == list(test_refl["intensity.prf.value"])
    assert list(new_rt["variance"]) == list(test_refl["intensity.prf.variance"])
    intstr = "sum"  # should apply partiality correction
    new_rt = choose_scaling_intensities(test_refl, intstr)
    assert list(new_rt["intensity"]) == list(
        test_refl["intensity.sum.value"] / test_refl["partiality"]
    )
    assert list(new_rt["variance"]) == pytest.approx(
        list(test_refl["intensity.sum.variance"] / (test_refl["partiality"] ** 2))
    )
    # If bad choice, currently return the prf values.
    intstr = "bad"
    new_rt = choose_scaling_intensities(test_refl, intstr)
    assert list(new_rt["intensity"]) == list(test_refl["intensity.prf.value"])
    assert list(new_rt["variance"]) == list(test_refl["intensity.prf.variance"])


def test_auto_scaling_model():
    params = generated_param()
    exp = generated_exp(scan=False)
    rt = generated_refl()
    params.__inject__("model", "auto")
    new_exp = create_auto_scaling_model(params, exp, [rt])
    assert new_exp[0].scaling_model.id_ == "KB"

    params = generated_param(absorption_term=True)
    exp = generated_exp(image_range=[1, 5])  # 5 degree wedge
    params.__inject__("model", "auto")
    new_exp = create_auto_scaling_model(params, exp, [rt])
    assert new_exp[0].scaling_model.id_ == "physical"
    assert len(new_exp[0].scaling_model.components["scale"].parameters) == 5
    assert len(new_exp[0].scaling_model.components["decay"].parameters) == 3
    assert "absorption" not in new_exp[0].scaling_model.components

    params = generated_param(absorption_term=True)
    exp = generated_exp(image_range=[1, 20])  # 20 degree wedge
    params.__inject__("model", "auto")
    new_exp = create_auto_scaling_model(params, exp, [rt])
    assert new_exp[0].scaling_model.id_ == "physical"
    assert len(new_exp[0].scaling_model.components["scale"].parameters) == 7
    assert len(new_exp[0].scaling_model.components["decay"].parameters) == 6
    assert "absorption" not in new_exp[0].scaling_model.components

    params = generated_param(absorption_term=True)
    exp = generated_exp(image_range=[1, 75])  # 20 degree wedge
    params.__inject__("model", "auto")
    new_exp = create_auto_scaling_model(params, exp, [rt])
    assert new_exp[0].scaling_model.id_ == "physical"
    assert len(new_exp[0].scaling_model.components["scale"].parameters) == 12
    assert len(new_exp[0].scaling_model.components["decay"].parameters) == 10
    assert "absorption" in new_exp[0].scaling_model.components

    # Now test overwrite_existing_models option
    params.overwrite_existing_models = True
    params.model = "KB"
    newer_exp = create_scaling_model(params, new_exp, [rt])
    assert isinstance(newer_exp[0].scaling_model, KBScalingModel)
