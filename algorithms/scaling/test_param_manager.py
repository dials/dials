"""
Tests for the active parameter manager module.
"""
from __future__ import absolute_import, division, print_function
import pytest
from mock import Mock
from scitbx import sparse
from dials.array_family import flex
from dials.algorithms.scaling.active_parameter_managers import (
    multi_active_parameter_manager,
    active_parameter_manager,
    ConcurrentAPMFactory,
    ConsecutiveAPMFactory,
)
from dials.algorithms.scaling.parameter_handler import (
    scaling_active_parameter_manager,
    create_apm_factory,
)


def mock_component():
    """Return a mock component of a general model."""
    component = Mock()
    component.parameters = flex.double([1.0])
    component.n_params = 1
    component.var_cov_matrix = sparse.matrix(1, 1)
    component.parameter_esds = None
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
    # dm.consecutive_refinement_order = None
    # dm.id_ = 'single'
    return dm


def mock_multiscaler(mock_data_managers):
    """Return a mock data manager of a general model."""
    multi_dm = Mock()
    multi_dm.active_scalers = mock_data_managers
    multi_dm.id_ = "multi"
    return multi_dm


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
    assert list(components["scale"].parameters) == [2.0]
    assert list(components["decay"].parameters) == [1.5]
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
    assert components["scale"].parameter_esds == flex.double([0.1])
    assert components["decay"].parameter_esds == flex.double([0.2])


def test_multi_apm():

    """Test for the general multi_active_parameter_manage class."""

    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    components_2 = {"scale": mock_component(), "decay": mock_component()}

    multi_apm = multi_active_parameter_manager(
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
    assert components_1["scale"].parameter_esds == flex.double([0.1])
    assert components_1["decay"].parameter_esds == flex.double([0.2])
    assert components_2["scale"].parameter_esds == flex.double([0.3])

    # Test setting var_cov matrices for each component.
    var_cov = flex.double([1.0, 0.5, 0.5, 0.5, 2.0, 0.5, 0.5, 0.5, 3.0])
    var_cov.reshape(flex.grid(3, 3))
    multi_apm.calculate_model_state_uncertainties(var_cov)
    assert components_1["scale"].var_cov_matrix[0, 0] == 1.0
    assert components_1["decay"].var_cov_matrix[0, 0] == 2.0
    assert components_2["scale"].var_cov_matrix[0, 0] == 3.0


def test_concurrent_apm_factory():
    """Test the apm factory for concurrent refinement."""
    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    data_manager = mock_data_manager(components_1)

    apm_factory = ConcurrentAPMFactory(
        [data_manager], apm_type=active_parameter_manager
    )
    apm = apm_factory.make_next_apm()
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" in apm.components_list

    # Test the multi-mode=False option
    apm_factory = ConcurrentAPMFactory(
        [data_manager], apm_type=active_parameter_manager, multi_mode=False
    )
    apm = apm_factory.make_next_apm()
    assert isinstance(apm, active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" in apm.components_list

    data_manager = mock_data_manager({})
    with pytest.raises(ValueError):
        apm_factory = ConcurrentAPMFactory(
            [data_manager], apm_type=active_parameter_manager
        )
    with pytest.raises(ValueError):
        apm_factory = ConcurrentAPMFactory(
            [data_manager], apm_type=active_parameter_manager, multi_mode=False
        )

    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }
    components_2 = {"1": mock_component(), "2": mock_component()}
    data_manager_1 = mock_data_manager(components_1)
    data_manager_2 = mock_data_manager(components_2)

    apm_factory = ConcurrentAPMFactory(
        [data_manager_1, data_manager_2], apm_type=active_parameter_manager
    )
    assert apm_factory.n_cycles == 1
    multi_apm = apm_factory.make_next_apm()
    assert isinstance(multi_apm, multi_active_parameter_manager)
    for apm in multi_apm.apm_list:
        assert isinstance(apm, active_parameter_manager)
    assert "scale" in multi_apm.apm_list[0].components_list
    assert "decay" in multi_apm.apm_list[0].components_list
    assert "absorption" in multi_apm.apm_list[0].components_list
    assert "1" in multi_apm.apm_list[1].components_list
    assert "2" in multi_apm.apm_list[1].components_list


def test_consecutive_apm_factory():
    """Test the apm factory for consecutive refinement."""
    components_1 = {
        "scale": mock_component(),
        "decay": mock_component(),
        "absorption": mock_component(),
    }

    data_manager = mock_data_manager(components_1)
    data_manager.consecutive_refinement_order = [["scale", "decay"], ["absorption"]]

    # Test single dataset case.
    apm_factory = ConsecutiveAPMFactory(
        [data_manager], apm_type=active_parameter_manager
    )
    assert apm_factory.n_cycles == 2
    apm = apm_factory.make_next_apm()
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" not in apm.components_list
    apm = apm_factory.make_next_apm()
    assert isinstance(apm, multi_active_parameter_manager)
    assert "scale" not in apm.components_list
    assert "decay" not in apm.components_list
    assert "absorption" in apm.components_list

    # Test the multi-mode=False option
    apm_factory = ConsecutiveAPMFactory(
        [data_manager], apm_type=active_parameter_manager, multi_mode=False
    )
    apm = apm_factory.make_next_apm()
    assert isinstance(apm, active_parameter_manager)
    assert "scale" in apm.components_list
    assert "decay" in apm.components_list
    assert "absorption" not in apm.components_list
    apm = apm_factory.make_next_apm()
    assert isinstance(apm, active_parameter_manager)
    assert "scale" not in apm.components_list
    assert "decay" not in apm.components_list
    assert "absorption" in apm.components_list

    # Test multi dataset case.
    components_2 = {"1": mock_component(), "2": mock_component()}
    data_manager_2 = mock_data_manager(components_2)
    data_manager_2.consecutive_refinement_order = [["1"], ["2"]]

    apm_factory = ConsecutiveAPMFactory(
        [data_manager, data_manager_2], apm_type=active_parameter_manager
    )
    assert apm_factory.n_cycles == 2
    multi_apm = apm_factory.make_next_apm()
    assert isinstance(multi_apm, multi_active_parameter_manager)
    apm_1 = multi_apm.apm_list[0]
    assert "scale" in apm_1.components_list
    assert "decay" in apm_1.components_list
    assert "absorption" not in apm_1.components_list
    assert multi_apm.apm_list[1].components_list == ["1"]
    multi_apm = apm_factory.make_next_apm()
    assert isinstance(multi_apm, multi_active_parameter_manager)
    assert multi_apm.apm_list[0].components_list == ["absorption"]
    assert multi_apm.apm_list[1].components_list == ["2"]

    # Test multi dataset case with different number of cycles for each data_manager.
    components_2 = {"1": mock_component()}
    data_manager_2 = mock_data_manager(components_2)
    data_manager_2.consecutive_refinement_order = [["1"]]
    apm_factory = ConsecutiveAPMFactory(
        [data_manager, data_manager_2], apm_type=active_parameter_manager
    )
    assert apm_factory.n_cycles == 2
    multi_apm = apm_factory.make_next_apm()
    assert isinstance(multi_apm, multi_active_parameter_manager)
    apm_1 = multi_apm.apm_list[0]
    assert "scale" in apm_1.components_list
    assert "decay" in apm_1.components_list
    assert "absorption" not in apm_1.components_list
    assert multi_apm.apm_list[1].components_list == ["1"]
    multi_apm = apm_factory.make_next_apm()
    assert isinstance(multi_apm, multi_active_parameter_manager)
    assert multi_apm.apm_list[0].components_list == ["absorption"]
    # Only change relative to previous test case.
    assert multi_apm.apm_list[1].components_list == []


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


def test_create_apm_factory():
    """Test that the create_apm_factory function correctly interprets the scaler.id
    and concurrent/consecutive option to give the correct APMFactory class."""

    # Concurrent single apm
    components_2 = {"1": mock_scaling_component(2), "2": mock_scaling_component(2)}
    scaler = mock_data_manager(components_2)
    scaler.id_ = "single"
    scaler.params.scaling_options.concurrent = True
    apm_factory = create_apm_factory(scaler)
    assert isinstance(apm_factory.apm, multi_active_parameter_manager)
    assert isinstance(apm_factory, ConcurrentAPMFactory)

    # Consecutive single apm
    scaler.params.scaling_options.concurrent = False
    scaler.consecutive_refinement_order = [["1"], ["2"]]
    apm_factory = create_apm_factory(scaler)
    # assert apm_factory.multi_mode is False
    assert isinstance(apm_factory, ConsecutiveAPMFactory)

    # Concurrent multi apm
    components_1 = {"1": mock_scaling_component(3), "2": mock_scaling_component(3)}
    components_2 = {"3": mock_scaling_component(2), "4": mock_scaling_component(2)}
    scaler1 = mock_data_manager(components_1)
    scaler1.consecutive_refinement_order = [["1"], ["2"]]
    scaler2 = mock_data_manager(components_2)
    scaler2.consecutive_refinement_order = [["3", "4"]]
    multiscaler = mock_multiscaler([scaler1, scaler2])

    for i in ["multi", "target"]:
        multiscaler.id_ = i
        multiscaler.params.scaling_options.concurrent = True
        apm_factory = create_apm_factory(multiscaler)
        assert apm_factory.multi_mode is True
        assert isinstance(apm_factory.apm, multi_active_parameter_manager)
        assert isinstance(apm_factory, ConcurrentAPMFactory)
        multiscaler.params.scaling_options.concurrent = False
        apm_factory = create_apm_factory(multiscaler)
        assert apm_factory.multi_mode is True
        assert isinstance(apm_factory, ConsecutiveAPMFactory)
    with pytest.raises(AssertionError):
        multiscaler.id_ = "1"
        apm_factory = create_apm_factory(multiscaler)
