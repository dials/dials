"""
Tests for the scaling model classes.
"""

from __future__ import annotations

import copy
from unittest.mock import MagicMock, Mock

import pytest

from libtbx import phil

from dials.algorithms.scaling.model.model import (
    ArrayScalingModel,
    DoseDecay,
    KBScalingModel,
    PhysicalScalingModel,
    ScalingModelBase,
    calc_n_param_from_bins,
    initialise_smooth_input,
)
from dials.array_family import flex
from dials.util.options import ArgumentParser


@pytest.fixture(scope="module")
def default_params():
    """Return the default parsed params phil scope."""
    return generated_param()


@pytest.fixture(scope="module")
def test_reflections():
    """Return a reflection table."""
    return generated_refl()


@pytest.fixture
def mock_exp():
    """Mock experiments object."""
    exp = Mock()
    exp.scan.get_oscillation.return_value = [0.0, 1.0]
    exp.scan.get_oscillation_range.return_value = [0, 90]
    exp.beam.get_s0.return_value = (0.0, 0.0, 1.01)
    exp.goniometer.get_rotation_axis.return_value = (1.0, 0.0, 0.0)
    return exp


def generated_refl():
    """Create a reflection table."""
    rt = flex.reflection_table()
    rt["xyzobs.px.value"] = flex.vec3_double([(0.1, 0.1, 0.1), (0.1, 0.1, 0.2)])
    rt["s1"] = flex.vec3_double([(0.1, 0.1, 0.1), (0.1, 0.1, 1.1)])
    rt["s1c"] = flex.vec3_double([(0.1, 0.1, 0.1), (0.1, 0.1, 1.1)])
    rt["s0c"] = flex.vec3_double([(0.0, 0.0, 1.0), (0.0, 0.0, 1.0)])
    rt["d"] = flex.double([1.0, 1.0])
    rt["batch"] = flex.int([0, 1])
    return rt


@pytest.fixture
def mock_errormodel():
    """Mock error model."""
    em = MagicMock()
    em.parameters = [1.0, 0.1]
    em.update_variances.return_value = flex.double([1.0, 1.1, 1.0, 1.0])
    return em


@pytest.fixture()
def mock_physical_params():
    """Return a mock params object for a physical model."""
    params = Mock()
    params.physical.scale_interval = 10.0
    params.physical.decay_correction = True
    params.physical.decay_interval = 15.0
    params.physical.absorption_correction = True
    params.physical.lmax = 4
    params.physical.decay_restraint = 1e-1
    params.physical.surface_weight = "auto"
    return params


def generated_param():
    """Generate the default scaling parameters object."""
    phil_scope = phil.parse(
        """
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.model.model.model_phil_scope
  """,
        process_includes=True,
    )

    parser = ArgumentParser(phil=phil_scope, check_format=False)
    parameters, _ = parser.parse_args(args=[], quick_parse=True, show_diff_phil=False)
    parameters.array.modulation_correction = True
    return parameters


def test_ScalingModelBase(mock_errormodel):
    """Test for base scaling model class"""

    configdict = {"corrections": ["1"]}
    SM_base = ScalingModelBase(configdict)
    assert not SM_base.is_scaled
    SM_base.set_scaling_model_as_scaled()
    assert SM_base.is_scaled
    SM_base.set_scaling_model_as_unscaled()
    assert not SM_base.is_scaled
    assert SM_base.configdict == configdict
    assert not SM_base.components
    _ = SM_base.to_dict()
    SM_base.set_error_model(mock_errormodel)
    assert SM_base.configdict["error_model_parameters"] == mock_errormodel.parameters
    assert SM_base.error_model is mock_errormodel
    assert "Scaling model" in str(SM_base)


def test_model_creation_from_data(default_params, mock_exp, test_reflections):
    """Test the factory creation of the three standard scaling models with the
    default params."""

    _ = KBScalingModel.from_data(default_params, [], [])

    _ = PhysicalScalingModel.from_data(default_params, mock_exp, test_reflections)

    _ = ArrayScalingModel.from_data(default_params, mock_exp, test_reflections)


def test_KBScalingModel():
    """Test for the KB Scaling Model."""

    # Test standard initialisation method.
    configdict = {"corrections": ["scale", "decay"]}
    parameters_dict = {
        "scale": {
            "parameters": flex.double([1.2]),
            "parameter_esds": flex.double([0.1]),
        },
        "decay": {
            "parameters": flex.double([0.01]),
            "parameter_esds": flex.double([0.02]),
        },
    }
    KBmodel = KBScalingModel(parameters_dict, configdict)
    assert KBmodel.id_ == "KB"
    assert "scale" in KBmodel.components
    assert "decay" in KBmodel.components
    assert list(KBmodel.components["scale"].parameters) == [1.2]
    assert list(KBmodel.components["decay"].parameters) == [0.01]
    assert list(KBmodel.components["scale"].parameter_esds) == [0.1]
    assert list(KBmodel.components["decay"].parameter_esds) == [0.02]

    # Test from_dict initialisation method.
    KB_dict = {
        "__id__": "KB",
        "scale": {
            "n_parameters": 1,
            "parameters": [0.5],
            "est_standard_devs": [0.05],
            "null_parameter_value": 1,
        },
        "configuration_parameters": {"corrections": ["scale"]},
    }
    KBmodel = KBScalingModel.from_dict(KB_dict)
    assert KBmodel.is_scaled is True
    assert "scale" in KBmodel.components
    assert "decay" not in KBmodel.components
    assert list(KBmodel.components["scale"].parameters) == [0.5]
    assert list(KBmodel.components["scale"].parameter_esds) == [0.05]

    new_dict = KBmodel.to_dict()
    assert new_dict == KB_dict

    # Test again with all parameters
    KB_dict = {
        "__id__": "KB",
        "scale": {
            "n_parameters": 1,
            "parameters": [0.5],
            "est_standard_devs": [0.05],
            "null_parameter_value": 1,
        },
        "decay": {
            "n_parameters": 1,
            "parameters": [0.2],
            "est_standard_devs": [0.02],
            "null_parameter_value": 0,
        },
        "configuration_parameters": {"corrections": ["scale", "decay"]},
    }
    KBmodel = KBScalingModel.from_dict(KB_dict)
    assert KBmodel.is_scaled is True
    assert "scale" in KBmodel.components
    assert "decay" in KBmodel.components
    assert list(KBmodel.components["scale"].parameters) == [0.5]
    assert list(KBmodel.components["scale"].parameter_esds) == [0.05]
    assert list(KBmodel.components["decay"].parameters) == [0.2]
    assert list(KBmodel.components["decay"].parameter_esds) == [0.02]

    new_dict = KBmodel.to_dict()
    assert new_dict == KB_dict

    with pytest.raises(RuntimeError):
        KB_dict["__id__"] = "physical"
        KBmodel = KBScalingModel.from_dict(KB_dict)

    assert KBmodel.consecutive_refinement_order == [["scale", "decay"]]
    assert "Decay component" in str(KBmodel)


def test_physical_model_from_data(mock_physical_params, mock_exp, test_reflections):
    """Test that it passes the correct dict to physical model."""
    physicalmodel = PhysicalScalingModel.from_data(
        mock_physical_params, mock_exp, test_reflections
    )
    assert physicalmodel.configdict["lmax"] == (mock_physical_params.physical.lmax)
    assert physicalmodel.components["absorption"].n_params == 24
    assert list(physicalmodel.components["absorption"].parameters) == [0.0] * 24

    # test updating the absorption parameters
    mock_physical_params.physical.absorption_level = "high"
    physicalmodel.update(mock_physical_params)
    assert len(physicalmodel.components["absorption"].parameters) == 48
    assert physicalmodel.configdict["abs_surface_weight"] == 5e3

    mock_physical_params.physical.absorption_level = "medium"
    physicalmodel.update(mock_physical_params)
    assert len(physicalmodel.components["absorption"].parameters) == 48
    assert physicalmodel.configdict["abs_surface_weight"] == 5e4

    mock_physical_params.physical.absorption_level = None
    mock_physical_params.physical.lmax = 4
    physicalmodel.update(mock_physical_params)
    assert len(physicalmodel.components["absorption"].parameters) == 24
    assert physicalmodel.configdict["abs_surface_weight"] == 5e4

    mock_physical_params.physical.surface_weight = 1e5
    physicalmodel.update(mock_physical_params)
    assert len(physicalmodel.components["absorption"].parameters) == 24
    assert physicalmodel.configdict["abs_surface_weight"] == 1e5

    # try fixing a parameter
    mock_physical_params.physical.correction.fix = ["decay"]
    physicalmodel = PhysicalScalingModel.from_data(
        mock_physical_params, mock_exp, test_reflections
    )
    assert physicalmodel.configdict["lmax"] == (mock_physical_params.physical.lmax)
    assert physicalmodel.components["absorption"].n_params == 24
    assert list(physicalmodel.components["absorption"].parameters) == [0.0] * 24
    assert physicalmodel.fixed_components == ["decay"]


def test_PhysicalScalingModel(test_reflections, mock_exp):
    """Test the PhysicalScalingModel class."""
    configdict = {
        "corrections": ["scale", "decay", "absorption"],
        "s_norm_fac": 1.0,
        "scale_rot_interval": 2.0,
        "d_norm_fac": 1.0,
        "decay_rot_interval": 2.0,
        "lmax": 1,
        "abs_surface_weight": 1e6,
    }

    parameters_dict = {
        "scale": {"parameters": flex.double([1.2, 1.1]), "parameter_esds": None},
        "decay": {"parameters": flex.double([0.1, 0.2]), "parameter_esds": None},
        "absorption": {
            "parameters": flex.double([0.01, 0.01, 0.01]),
            "parameter_esds": None,
        },
    }

    # Test standard factory initialisation
    physicalmodel = PhysicalScalingModel(parameters_dict, configdict)
    assert physicalmodel.id_ == "physical"
    assert "Absorption component" in str(physicalmodel)
    comps = physicalmodel.components
    assert "scale" in comps
    assert "absorption" in comps
    assert "decay" in comps
    assert list(comps["scale"].parameters) == [1.2, 1.1]
    assert list(comps["decay"].parameters) == [0.1, 0.2]
    assert list(comps["absorption"].parameters) == [0.01, 0.01, 0.01]

    # Test configure reflection table
    mock_params = Mock()
    mock_params.physical.decay_restraint = 0.0
    physicalmodel.configure_components(test_reflections, mock_exp, mock_params)

    # Test from_dict initialisation method.
    physical_dict = {
        "__id__": "physical",
        "scale": {
            "n_parameters": 2,
            "parameters": [0.5, 1.0],
            "est_standard_devs": [0.05, 0.1],
            "null_parameter_value": 1,
        },
        "configuration_parameters": {
            "corrections": ["scale"],
            "s_norm_fac": 0.1,
            "scale_rot_interval": 10.0,
            "decay_restaint": 1e-1,
        },
    }
    physicalmodel = PhysicalScalingModel.from_dict(physical_dict)
    assert physicalmodel.id_ == "physical"
    assert "scale" in physicalmodel.components
    assert "absorption" not in physicalmodel.components
    assert "decay" not in physicalmodel.components
    assert list(physicalmodel.components["scale"].parameters) == [0.5, 1.0]
    assert list(physicalmodel.components["scale"].parameter_esds) == [0.05, 0.1]

    new_dict = physicalmodel.to_dict()
    assert new_dict == physical_dict

    # Test from_dict initialisation method for all components.
    physical_dict = {
        "__id__": "physical",
        "scale": {
            "n_parameters": 2,
            "parameters": [0.5, 1.0],
            "est_standard_devs": [0.05, 0.1],
            "null_parameter_value": 1,
        },
        "decay": {
            "n_parameters": 2,
            "parameters": [0.1, 0.2],
            "est_standard_devs": [0.01, 0.01],
            "null_parameter_value": 0,
        },
        "absorption": {
            "n_parameters": 3,
            "parameters": [0.0, 0.1, 0.2],
            "est_standard_devs": [0.01, 0.02, 0.03],
            "null_parameter_value": 0,
        },
        "configuration_parameters": {
            "corrections": ["scale", "decay", "absorption"],
            "s_norm_fac": 0.1,
            "scale_rot_interval": 10.0,
            "d_norm_fac": 0.2,
            "decay_rot_interval": 20.0,
            "lmax": 1,
            "abs_surface_weight": 1e6,
        },
    }
    physicalmodel = PhysicalScalingModel.from_dict(physical_dict)
    assert physicalmodel.id_ == "physical"
    assert "scale" in physicalmodel.components
    assert "absorption" in physicalmodel.components
    assert "decay" in physicalmodel.components
    assert list(physicalmodel.components["scale"].parameters) == [0.5, 1.0]
    assert list(physicalmodel.components["scale"].parameter_esds) == [0.05, 0.1]
    assert list(physicalmodel.components["decay"].parameters) == [0.1, 0.2]
    assert list(physicalmodel.components["decay"].parameter_esds) == [0.01, 0.01]
    assert list(physicalmodel.components["absorption"].parameters) == [0.0, 0.1, 0.2]
    assert list(physicalmodel.components["absorption"].parameter_esds) == [
        0.01,
        0.02,
        0.03,
    ]

    new_dict = physicalmodel.to_dict()
    assert new_dict == physical_dict

    with pytest.raises(RuntimeError):
        physical_dict["__id__"] = "array"
        physicalmodel = PhysicalScalingModel.from_dict(physical_dict)

    assert len(physicalmodel.consecutive_refinement_order) == 2
    assert "Absorption component:" in str(physicalmodel)

    # test limit batch range
    parameters_dict = {
        "scale": {
            "n_parameters": 11,
            "parameters": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5],
            "parameter_esds": None,
        },
        "decay": {
            "n_parameters": 11,
            "parameters": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1],
            "parameter_esds": None,
        },
    }
    configdict = {
        "corrections": ["scale", "decay"],
        "s_norm_fac": 0.1,
        "scale_rot_interval": 10.0,
        "d_norm_fac": 0.1,
        "decay_rot_interval": 10.0,
        "valid_image_range": (1, 90),
        "valid_osc_range": (0, 90),
    }
    physical = PhysicalScalingModel(parameters_dict, configdict)
    physical.limit_image_range((1, 50))
    assert list(physical.components["scale"].parameters) == [
        0.5,
        1.0,
        1.5,
        2.0,
        2.5,
        3.0,
        3.5,
    ]
    assert list(physical.components["decay"].parameters) == [
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
    ]
    assert physical.configdict["valid_osc_range"] == (0, 50)
    assert physical.configdict["valid_image_range"] == (1, 50)

    # try edge cases
    # if restricted by > rot int, then reduce number of params and shift offset
    # if necessary
    physical = PhysicalScalingModel(
        copy.deepcopy(parameters_dict), copy.deepcopy(configdict)
    )
    physical.limit_image_range((7, 45))
    assert list(physical.components["scale"].parameters) == [
        1.0,
        1.5,
        2.0,
        2.5,
        3.0,
        3.5,
    ]
    assert list(physical.components["decay"].parameters) == [
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
    ]
    assert physical.configdict["valid_osc_range"] == (6, 45)
    assert physical.configdict["valid_image_range"] == (7, 45)

    # if not restricted by > rot int, then should 'squeeze' the parameters closer
    # in rotation angle, leaving the same number of params (as reducing number of
    # params would give parameter spacing greater than initially specified rot int)
    physical = PhysicalScalingModel(
        copy.deepcopy(parameters_dict), copy.deepcopy(configdict)
    )
    physical.limit_image_range((5, 45))
    assert list(physical.components["scale"].parameters) == [
        0.5,
        1.0,
        1.5,
        2.0,
        2.5,
        3.0,
        3.5,
    ]
    assert list(physical.components["decay"].parameters) == [
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
    ]
    assert physical.configdict["valid_osc_range"] == (4, 45)
    assert physical.configdict["valid_image_range"] == (5, 45)


def test_DoseDecayModel(test_reflections, mock_exp):
    """Test the DoseDecay class."""
    configdict = {
        "corrections": ["scale", "decay", "relative_B"],
        "s_norm_fac": 1.0,
        "scale_rot_interval": 2.0,
        "resolution_dependence": "quadratic",
    }

    parameters_dict = {
        "scale": {"parameters": flex.double([1.2, 1.1]), "parameter_esds": None},
        "decay": {"parameters": flex.double([0.1]), "parameter_esds": None},
        "relative_B": {"parameters": flex.double([-0.1]), "parameter_esds": None},
    }

    def test_setup(model):
        assert model.id_ == "dose_decay"
        assert set(model.components.keys()) == {"scale", "decay", "relative_B"}
        assert list(model.components["scale"].parameters) == [1.2, 1.1]
        assert list(model.components["decay"].parameters) == [0.1]
        assert list(model.components["relative_B"].parameters) == [-0.1]

    # Test standard factory initialisation
    model = DoseDecay(parameters_dict, configdict)
    test_setup(model)

    # Test configure reflection table - adds data to components
    mock_params = Mock()
    model.configure_components(test_reflections, mock_exp, mock_params)
    assert list(model.components["decay"].data["d"]) == list(test_reflections["d"])

    d = model.to_dict()
    model_from_dict = DoseDecay.from_dict(d)
    test_setup(model_from_dict)


def test_ArrayScalingModel(test_reflections, mock_exp):
    """Test the ArrayScalingModel class."""

    configdict = {
        "corrections": ["decay", "absorption", "modulation"],
        "n_res_param": 2,
        "n_time_param": 2,
        "resmin": 1.0,
        "res_bin_width": 1.0,
        "time_norm_fac": 1.0,
        "time_rot_interval": 1.0,
        "n_x_param": 2,
        "n_y_param": 2,
        "xmin": 0.0,
        "ymin": 0.0,
        "x_bin_width": 1.0,
        "y_bin_width": 2.0,
        "n_x_mod_param": 2,
        "n_y_mod_param": 2,
        "x_det_bin_width": 2.0,
        "y_det_bin_width": 2.0,
    }

    parameters_dict = {
        "decay": {
            "parameters": flex.double([1.2, 1.1, 1.0, 0.9]),
            "parameter_esds": None,
        },
        "absorption": {
            "parameters": flex.double([0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2]),
            "parameter_esds": None,
        },
        "modulation": {"parameters": flex.double(4, 0.9), "parameter_esds": None},
    }

    # Test standard factory initialisation
    arraymodel = ArrayScalingModel(parameters_dict, configdict)
    assert arraymodel.id_ == "array"
    assert "decay" in arraymodel.components
    assert "absorption" in arraymodel.components
    assert "modulation" in arraymodel.components
    assert list(arraymodel.components["decay"].parameters) == [1.2, 1.1, 1.0, 0.9]
    assert list(arraymodel.components["absorption"].parameters) == [
        0.1,
        0.2,
        0.1,
        0.2,
        0.1,
        0.2,
        0.1,
        0.2,
    ]
    assert list(arraymodel.components["modulation"].parameters) == 4 * [0.9]

    # Test configure reflection table
    _ = arraymodel.configure_components(test_reflections, mock_exp, [])

    # Test from_dict initialisation method for previous model case.
    init_dict = arraymodel.to_dict()
    new_array_model = ArrayScalingModel.from_dict(init_dict)
    assert new_array_model.id_ == "array"
    comps = new_array_model.components
    assert "modulation" in comps
    assert "absorption" in comps
    assert "decay" in comps
    assert list(comps["decay"].parameters) == [1.2, 1.1, 1.0, 0.9]
    assert list(comps["absorption"].parameters) == [
        0.1,
        0.2,
        0.1,
        0.2,
        0.1,
        0.2,
        0.1,
        0.2,
    ]
    assert list(comps["modulation"].parameters) == 4 * [0.9]

    # Test from_dict initialisation method for another case.
    array_dict = {
        "__id__": "array",
        "decay": {
            "n_parameters": 4,
            "parameters": [0.5, 1.0, 0.4, 1.0],
            "null_parameter_value": 1.0,
            "est_standard_devs": [0.05, 0.1, 0.05, 0.1],
        },
        "configuration_parameters": {
            "corrections": ["decay"],
            "n_res_param": 2,
            "n_time_param": 2,
            "resmin": 1.0,
            "res_bin_width": 1.0,
            "time_norm_fac": 1.0,
            "time_rot_interval": 1.0,
        },
    }
    arraymodel = ArrayScalingModel.from_dict(array_dict)
    assert arraymodel.id_ == "array"
    comps = arraymodel.components
    assert "modulation" not in comps
    assert "absorption" not in comps
    assert "decay" in comps
    assert list(comps["decay"].parameters) == [0.5, 1.0, 0.4, 1.0]
    assert list(comps["decay"].parameter_esds) == [0.05, 0.1, 0.05, 0.1]

    new_dict = arraymodel.to_dict()
    assert new_dict == array_dict

    with pytest.raises(RuntimeError):
        array_dict["__id__"] = "physical"
        arraymodel = ArrayScalingModel.from_dict(array_dict)

    assert len(arraymodel.consecutive_refinement_order) == 3
    assert "Decay component" in str(arraymodel)

    # test limit batch range
    configdict = {
        "corrections": ["decay", "absorption"],
        "n_res_param": 2,
        "n_time_param": 3,
        "resmin": 1.0,
        "res_bin_width": 1.0,
        "time_norm_fac": 0.1,
        "time_rot_interval": 10.0,
        "n_x_param": 2,
        "n_y_param": 2,
        "xmin": 0.0,
        "ymin": 0.0,
        "x_bin_width": 1.0,
        "y_bin_width": 2.0,
        "n_x_mod_param": 2,
        "n_y_mod_param": 2,
        "x_det_bin_width": 2.0,
        "y_det_bin_width": 2.0,
        "valid_image_range": (1, 20),
        "valid_osc_range": (0, 20),
    }

    parameters_dict = {
        "decay": {
            "parameters": flex.double([1.2, 1.1, 1.0, 0.9, 0.8, 0.7]),
            "parameter_esds": None,
        },
        "absorption": {
            "parameters": flex.double(
                [0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.3, 0.4, 0.3, 0.4]
            ),
            "parameter_esds": None,
        },
    }
    array = ArrayScalingModel(parameters_dict, configdict)
    array.limit_image_range((1, 10))
    assert list(array.components["decay"].parameters) == [1.2, 1.1, 1.0, 0.9]
    assert list(array.components["absorption"].parameters) == [
        0.1,
        0.2,
        0.1,
        0.2,
        0.1,
        0.2,
        0.1,
        0.2,
    ]
    assert array.configdict["n_time_param"] == 2
    assert array.configdict["valid_image_range"] == (1, 10)
    assert array.configdict["valid_osc_range"] == (0, 10)


def test_model_factory_utilities():
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
    # rotation interval in degrees, returning
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
    assert (n_param, norm_fac, rot_int) == (5, 0.999 * 3.0 / 10.0, 10.0 / 3.0)
    n_param, norm_fac, rot_int = initialise_smooth_input([0, 10], 2.0, 4.99)
    assert (n_param, norm_fac, rot_int) == (5, 0.999 * 6.0 / 10.0, 10.0 / 3.0)
