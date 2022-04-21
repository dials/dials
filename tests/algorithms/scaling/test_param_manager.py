"""
Tests for the active parameter manager module.
"""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from scitbx import sparse

from dials.algorithms.scaling.active_parameter_managers import (
    ParameterManagerGenerator,
    active_parameter_manager,
    multi_active_parameter_manager,
    shared_active_parameter_manager,
)
from dials.algorithms.scaling.parameter_handler import (
    ScalingParameterManagerGenerator,
    scaling_active_parameter_manager,
)
from dials.algorithms.scaling.target_function import ScalingTarget
from dials.array_family import flex


def mock_component():
    """Return a mock component of a general model."""
    component = Mock()
    component.free_parameters = flex.double([1.0])
    component.free_parameter_esds = None
    component.n_params = 1
    component.var_cov_matrix = sparse.matrix(1, 1)
    return component


def mock_scaling_component(n_refl):
    """Return a mock component of a general model."""
    component = mock_component()
    component.calculate_scales.return_value = flex.double(n_refl, 1.0)
    component.n_refl = [n_refl]
    return component


def mock_data_manager(components):
    """Return a mock data manager of a general model."""
    dm = Mock()
    dm.components = components
    dm.fixed_components = []
    return dm


def test_general_apm():
    """Test for the general active_parameter_manage class."""
    components = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }

    apm = active_parameter_manager(components, ["scale", "decay"])
    assert "decay" in apm.components_list
    assert "scale" in apm.components_list
    assert "absorption" not in apm.components_list
    assert apm.n_active_params == (
        components["scale"].n_params + components["decay"].n_params
    )
    n_cumul = 0
    for component in apm.components:
        assert apm.components[component]["n_params"] == components[component].n_params
        assert apm.components[component]["start_idx"] == n_cumul
        assert (
            apm.components[component]["end_idx"]
            == n_cumul + apm.components[component]["n_params"]
        )
        n_cumul += apm.components[component]["n_params"]

    apm.set_param_vals(flex.double([2.0, 1.5]))
    assert apm.get_param_vals() == flex.double([2.0, 1.5])
    # Test params were updated in components
    assert list(components["scale"].free_parameters) == [2.0]
    assert list(components["decay"].free_parameters) == [1.5]
    # Test selection of parameters
    decay_params = apm.select_parameters("decay")
    assert len(decay_params) == 1
    assert decay_params[0] == 1.5

    # Test calculate model state uncertainties
    var_cov = flex.double([1.0, 0.5, 0.5, 2.0])
    var_cov.reshape(flex.grid(2, 2))
    apm.calculate_model_state_uncertainties(var_cov)
    assert components["scale"].var_cov_matrix[0, 0] == 1.0
    assert components["decay"].var_cov_matrix[0, 0] == 2.0

    # Test set param esds.
    apm.set_param_esds(flex.double([0.1, 0.2]))
    assert components["scale"].free_parameter_esds == flex.double([0.1])
    assert components["decay"].free_parameter_esds == flex.double([0.2])


def test_multi_apm():
    """Test for the general multi_active_parameter_manage class."""

    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    components_2 = {"scale": mock_component(), "decay": mock_component()}

    multi_apm = multi_active_parameter_manager(
        ScalingTarget(),
        [components_1, components_2],
        [["scale", "decay"], ["scale"]],
        active_parameter_manager,
    )

    # Test correct setup of apm_list attribute.
    for apm in multi_apm.apm_list:
        assert isinstance(apm, active_parameter_manager)
    assert len(multi_apm.apm_list) == 2
    assert multi_apm.components_list == ["scale", "decay", "scale"]
    assert multi_apm.n_active_params == 3
    assert multi_apm.apm_data[0] == {"start_idx": 0, "end_idx": 2}
    assert multi_apm.apm_data[1] == {"start_idx": 2, "end_idx": 3}

    # Test parameter selection.
    multi_apm.set_param_vals(flex.double([3.0, 2.5, 2.0]))
    assert multi_apm.get_param_vals() == flex.double([3.0, 2.5, 2.0])
    assert multi_apm.select_parameters(0) == flex.double([3.0, 2.5])
    assert multi_apm.select_parameters(1) == flex.double([2.0])

    # Test setting parameter esds.
    multi_apm.set_param_esds(flex.double([0.1, 0.2, 0.3]))
    assert components_1["scale"].free_parameter_esds == flex.double([0.1])
    assert components_1["decay"].free_parameter_esds == flex.double([0.2])
    assert components_2["scale"].free_parameter_esds == flex.double([0.3])

    # Test setting var_cov matrices for each component.
    var_cov = flex.double([1.0, 0.5, 0.5, 0.5, 2.0, 0.5, 0.5, 0.5, 3.0])
    var_cov.reshape(flex.grid(3, 3))
    multi_apm.calculate_model_state_uncertainties(var_cov)
    assert components_1["scale"].var_cov_matrix[0, 0] == 1.0
    assert components_1["decay"].var_cov_matrix[0, 0] == 2.0
    assert components_2["scale"].var_cov_matrix[0, 0] == 3.0


def test_shared_apm():
    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    components_2 = {"scale": mock_component(), "decay": mock_component()}

    multi_apm = shared_active_parameter_manager(
        ScalingTarget(),
        [components_1, components_2],
        [["scale", "decay"], ["scale", "decay"]],
        active_parameter_manager,
        shared="decay",
    )

    # Test correct setup of apm_list attribute.
    for apm in multi_apm.apm_list:
        assert isinstance(apm, active_parameter_manager)
    assert len(multi_apm.apm_list) == 2
    assert multi_apm.components_list == ["scale", "decay", "scale", "decay"]
    assert multi_apm.n_active_params == 4
    # assert multi_apm.apm_data[0] == {"start_idx": 0, "end_idx": 2}
    # assert multi_apm.apm_data[1] == {"start_idx": 2, "end_idx": 4}
    assert list(multi_apm.reducing_matrix.as_dense_matrix()) == [
        1,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        1,
        0,
        1,
        0,
    ]
    assert multi_apm.apm_data[0]["start_idx"] == 0
    assert multi_apm.apm_data[0]["end_idx"] == 2
    assert multi_apm.apm_data[1]["start_idx"] == 2
    assert multi_apm.apm_data[1]["end_idx"] == 4
    assert list(multi_apm.apm_data[0]["apm_sel"]) == [0, 1]
    assert list(multi_apm.apm_data[1]["apm_sel"]) == [2, 1]


def test_ParameterManagerGenerator_concurrent():
    """Test the apm factory for concurrent refinement."""
    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    data_manager = mock_data_manager(components_1)

    pmg = ParameterManagerGenerator(
        [data_manager],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="concurrent",
    )
    apms = pmg.parameter_managers()
    assert len(apms) == 1
    apm = apms[0]
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" in apm.components_list

    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    components_2 = {"1": mock_component(), "2": mock_component()}
    data_manager_1 = mock_data_manager(components_1)
    data_manager_2 = mock_data_manager(components_2)

    pmg = ParameterManagerGenerator(
        [data_manager_1, data_manager_2],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="concurrent",
    )
    multi_apms = pmg.parameter_managers()
    assert len(multi_apms) == 1
    multi_apm = multi_apms[0]
    assert isinstance(multi_apm, multi_active_parameter_manager)
    for apm in multi_apm.apm_list:
        assert isinstance(apm, active_parameter_manager)
    assert "scale" in multi_apm.apm_list[0].components_list
    assert "decay" in multi_apm.apm_list[0].components_list
    assert "absorption" in multi_apm.apm_list[0].components_list
    assert "1" in multi_apm.apm_list[1].components_list
    assert "2" in multi_apm.apm_list[1].components_list

    # now try fixing a component
    data_manager.fixed_components = ["absorption"]
    pmg = ParameterManagerGenerator(
        [data_manager],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="concurrent",
    )
    apms = pmg.parameter_managers()
    assert len(apms) == 1
    apm = apms[0]
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" not in apm.components_list


def test_ParameterManagerGenerator_consecutive():
    """Test the apm factory for consecutive refinement."""
    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }

    data_manager = mock_data_manager(components_1)
    data_manager.consecutive_refinement_order = [["scale", "decay"], ["absorption"]]

    # Test single dataset case.
    pmg = ParameterManagerGenerator(
        [data_manager],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="consecutive",
    )
    apms = list(pmg.parameter_managers())
    assert len(apms) == 2
    apm = apms[0]
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" not in apm.components_list
    apm = apms[1]
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" not in apm.components_list
    assert "decay" not in apm.components_list
    assert "absorption" in apm.components_list

    # Test multi dataset case.
    components_2 = {"1": mock_component(), "2": mock_component()}
    data_manager_2 = mock_data_manager(components_2)
    data_manager_2.consecutive_refinement_order = [["1"], ["2"]]

    pmg = ParameterManagerGenerator(
        [data_manager, data_manager_2],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="consecutive",
    )
    apms = list(pmg.parameter_managers())
    assert len(apms) == 2
    multi_apm = apms[0]
    assert isinstance(multi_apm, multi_active_parameter_manager)
    apm_1 = multi_apm.apm_list[0]
    assert "scale" in apm_1.components_list
    assert "decay" in apm_1.components_list
    assert "absorption" not in apm_1.components_list
    assert multi_apm.apm_list[1].components_list == ["1"]
    multi_apm = apms[1]
    assert isinstance(multi_apm, multi_active_parameter_manager)
    assert multi_apm.apm_list[0].components_list == ["absorption"]
    assert multi_apm.apm_list[1].components_list == ["2"]

    # Test multi dataset case with different number of cycles for each data_manager.
    components_2 = {"1": mock_component()}
    data_manager_2 = mock_data_manager(components_2)
    data_manager_2.consecutive_refinement_order = [["1"], ["2"]]
    pmg = ParameterManagerGenerator(
        [data_manager, data_manager_2],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="consecutive",
    )
    assert pmg.param_lists[0] == [["scale", "decay"], ["absorption"]]
    assert pmg.param_lists[1] == [["1"]]
    apms = list(pmg.parameter_managers())
    assert len(apms) == 2
    multi_apm = apms[0]
    assert isinstance(multi_apm, multi_active_parameter_manager)
    apm_1 = multi_apm.apm_list[0]
    assert "scale" in apm_1.components_list
    assert "decay" in apm_1.components_list
    assert "absorption" not in apm_1.components_list
    assert multi_apm.apm_list[1].components_list == ["1"]
    multi_apm = apms[1]
    assert isinstance(multi_apm, multi_active_parameter_manager)
    assert multi_apm.apm_list[0].components_list == ["absorption"]
    # Only change relative to previous test case.
    assert multi_apm.apm_list[1].components_list == []

    # Test fixing the decay parameter.
    data_manager.fixed_components = ["decay"]
    pmg = ParameterManagerGenerator(
        [data_manager],
        apm_type=active_parameter_manager,
        target=ScalingTarget(),
        mode="consecutive",
    )
    apms = list(pmg.parameter_managers())
    assert len(apms) == 2
    apm = apms[0]
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" not in apm.components_list
    assert "absorption" not in apm.components_list
    apm = apms[1]
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" not in apm.components_list
    assert "decay" not in apm.components_list
    assert "absorption" in apm.components_list


def test_scaling_active_parameter_manager():
    """Test the scaling-specific parameter manager."""
    components_2 = {"1": mock_scaling_component(2), "2": mock_scaling_component(2)}
    scaling_apm = scaling_active_parameter_manager(components_2, ["1"])
    assert list(scaling_apm.constant_g_values[0]) == list(
        components_2["2"].calculate_scales()
    )
    assert len(scaling_apm.constant_g_values) == 1
    assert scaling_apm.n_obs == [2]

    # Test that no constant_g_values if both components selected
    scaling_apm = scaling_active_parameter_manager(components_2, ["1", "2"])
    assert scaling_apm.constant_g_values is None

    # Check that one can't initialise with an unequal number of reflections,
    # either within the selection or overall.
    with pytest.raises(AssertionError):
        components_2 = {"1": mock_scaling_component(2), "2": mock_scaling_component(1)}
        scaling_apm = scaling_active_parameter_manager(components_2, ["1", "2"])
    with pytest.raises(AssertionError):
        components_2 = {"1": mock_scaling_component(2), "2": mock_scaling_component(1)}
        scaling_apm = scaling_active_parameter_manager(components_2, ["1"])

    data_manager = mock_data_manager(components_2)
    pmg = ScalingParameterManagerGenerator(
        [data_manager], target=ScalingTarget(), mode="concurrent"
    )
    assert isinstance(pmg.apm_type, type(scaling_active_parameter_manager))
