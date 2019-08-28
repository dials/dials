"""
Tests for the scaling model classes.
"""
from __future__ import absolute_import, division, print_function
import copy
import pytest
from mock import Mock, MagicMock
from dials.array_family import flex
from dials.algorithms.scaling.model.model import (
    ScalingModelBase,
    KBScalingModel,
    PhysicalScalingModel,
    ArrayScalingModel,
)  # , map_old_to_new_range


@pytest.fixture(scope="module")
def test_reflections():
    """Return a reflection table."""
    return generated_refl()


@pytest.fixture
def mock_exp():
    """Mock experiments object."""
    exp = Mock()
    exp.scan.get_oscillation.return_value = [0.0, 1.0]
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
    em.refined_parameters = [1.0, 0.1]
    em.update_variances.return_value = flex.double([1.0, 1.1, 1.0, 1.0])
    return em


def test_ScalingModelBase(mock_errormodel):
    """Test for base scaling model class"""

    class SM_base_filler(ScalingModelBase):
        """Fill in abstract methid"""

        def consecutive_refinement_order(self):
            """Fill in abstract method."""

        def from_dict(self, obj):
            """Fill in abstract method."""

        def configure_components(self, _, __, ___):
            """Fill in abstract method."""

    SM_base = SM_base_filler(configdict={})
    assert not SM_base.is_scaled
    SM_base.set_scaling_model_as_scaled()
    assert SM_base.is_scaled
    SM_base.set_scaling_model_as_unscaled()
    assert not SM_base.is_scaled
    assert SM_base.configdict == {}
    assert not SM_base.components
    _ = SM_base.to_dict()
    SM_base.set_error_model(mock_errormodel)
    assert (
        SM_base.configdict["error_model_parameters"]
        == mock_errormodel.refined_parameters
    )
    assert SM_base.error_model is mock_errormodel
    _ = SM_base.configure_components(1.0, 2.0, 0.0)  # Check method exists
    SM_base.show()


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
        "is_scaled": True,
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
        "is_scaled": True,
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
    KBmodel.show()


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
    physicalmodel.show()  # test show works with no parameter esds.
    comps = physicalmodel.components
    assert "scale" in comps
    assert "absorption" in comps
    assert "decay" in comps
    assert list(comps["scale"].parameters) == [1.2, 1.1]
    assert list(comps["decay"].parameters) == [0.1, 0.2]
    assert list(comps["absorption"].parameters) == [0.01, 0.01, 0.01]

    # Test configure reflection table
    mock_params = Mock()
    mock_params.parameterisation.decay_restraint = 0.0
    physicalmodel.configure_components(test_reflections, mock_exp, mock_params)
    # Test normalise components.
    physicalmodel.components["scale"].update_reflection_data()
    physicalmodel.components["scale"].calculate_scales_and_derivatives()
    physicalmodel.components["decay"].update_reflection_data()
    physicalmodel.components["decay"].calculate_scales_and_derivatives()
    physicalmodel.normalise_components()
    assert list(physicalmodel.components["scale"].parameters) == pytest.approx(
        [1.0091065, 0.925014], 1e-4
    )

    # Test from_dict initialisation method.
    physical_dict = {
        "__id__": "physical",
        "is_scaled": True,
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
        "is_scaled": True,
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
    physicalmodel.show()

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
        "is_scaled": True,
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
    arraymodel.show()

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
