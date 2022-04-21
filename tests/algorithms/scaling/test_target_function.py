"""
Test for the target function module.
"""

from __future__ import annotations

import copy
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pytest
from scipy.sparse import csc_matrix

from dxtbx.model import Crystal, Experiment, Scan
from libtbx import phil
from scitbx import sparse

from dials.algorithms.scaling.active_parameter_managers import (
    multi_active_parameter_manager,
)
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.model.model import PhysicalScalingModel
from dials.algorithms.scaling.parameter_handler import scaling_active_parameter_manager
from dials.algorithms.scaling.scaler import SingleScaler
from dials.algorithms.scaling.target_function import ScalingTarget, ScalingTargetFixedIH
from dials.array_family import flex
from dials.util.options import ArgumentParser


@pytest.fixture
def large_reflection_table():
    """Generate reflection table to test the basis and target function."""
    # these miller_idx/d_values don't make physical sense, but I didn't want to
    # have to write the tests for lots of reflections.
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double(
        [75.0, 10.0, 100.0, 25.0, 50.0, 100.0, 25.0, 20.0, 300.0, 10.0]
    )
    reflections["variance"] = flex.double(
        [50.0, 10.0, 100.0, 50.0, 10.0, 100.0, 50.0, 10.0, 100.0, 10.0]
    )
    reflections["inverse_scale_factor"] = flex.double(10, 1.0)
    reflections["miller_index"] = flex.miller_index(
        [
            (1, 0, 0),
            (0, 0, 1),
            (1, 0, 0),
            (1, 0, 0),
            (0, 0, 1),
            (1, 0, 0),
            (0, 4, 0),
            (0, 0, 1),
            (1, 0, 0),
            (0, 4, 0),
        ]
    )  # don't change
    reflections["d"] = flex.double(
        [2.0, 0.8, 2.0, 2.0, 0.8, 2.0, 2.0, 0.8, 2.0, 1.0]
    )  # don't change
    reflections["partiality"] = flex.double(10, 1.0)
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [(0.0, 0.0, i) for i in range(0, 45, 5)] + [(0.0, 0.0, 59.0)]
    )
    reflections["s1"] = flex.vec3_double([(0.0, 0.1, 1.0)] * 2 + [(0.0, 0.1, 20.0)] * 8)
    return reflections


@pytest.fixture
def small_reflection_table():
    """Generate reflection table to test the basis and target function."""
    # these miller_idx/d_values don't make physical sense, but I didn't want to
    # have to write the tests for lots of reflections.
    reflections = flex.reflection_table()
    reflections["intensity"] = flex.double([75.0, 10.0, 100.0, 50.0, 65.0])
    reflections["variance"] = flex.double([50.0, 10.0, 100.0, 50.0, 65.0])
    reflections["inverse_scale_factor"] = flex.double(5, 1.0)
    reflections["miller_index"] = flex.miller_index(
        [(1, 0, 0), (0, 0, 1), (1, 0, 0), (0, 0, 1), (0, 0, 2)]
    )  # don't change
    reflections["d"] = flex.double([2.0, 0.8, 2.0, 0.8, 1.5])  # don't change
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [(0.0, 0.0, i) for i in [0.0, 5.0, 10.0, 2.0, 8.0]]
    )
    reflections["s1"] = flex.vec3_double([(0.0, 0.1, 1.0)] * 2 + [(0.0, 0.1, 20.0)] * 3)
    return reflections


@pytest.fixture
def single_exp():
    """Generate an experiment object."""
    crystal = Crystal.from_dict(
        {
            "__id__": "crystal",
            "real_space_a": [1.0, 0.0, 0.0],
            "real_space_b": [0.0, 1.0, 0.0],
            "real_space_c": [0.0, 0.0, 2.0],
            "space_group_hall_symbol": " C 2y",
        }
    )
    scan = Scan(image_range=[0, 60], oscillation=[0.0, 1.0])
    return Experiment(scan=scan, crystal=crystal)


@pytest.fixture
def physical_param():
    """Generate the scaling phil param scope."""
    phil_scope = phil.parse(
        """
      include scope dials.algorithms.scaling.model.model.model_phil_scope
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  """,
        process_includes=True,
    )

    parser = ArgumentParser(phil=phil_scope, check_format=False)
    parameters, _ = parser.parse_args(args=[], quick_parse=True, show_diff_phil=False)
    parameters.model = "physical"
    parameters.physical.absorption_correction = False
    return parameters


def mock_single_Ih_table():
    """Mock Ih table to use for testing the target function."""
    Ih_table = Mock()
    Ih_table.inverse_scale_factors = np.array([1.0, 1.0 / 1.1, 1.0])
    Ih_table.intensities = np.array([10.0, 10.0, 12.0])
    Ih_table.Ih_values = np.array([11.0, 11.0, 11.0])
    # These values should give residuals of [-1.0, 0.0, 1.0]
    Ih_table.weights = np.array([1.0, 1.0, 1.0])
    Ih_table.size = 3
    Ih_table.derivatives = sparse.matrix(3, 1, [{0: 1.0, 1: 2.0, 2: 3.0}])
    Ih_table.h_index_matrix = sparse.matrix(3, 2, [{0: 1, 1: 1}, {2: 1}])
    Ih_table.h_expand_matrix = Ih_table.h_index_matrix.transpose()
    Ih_table._csc_h_index_matrix = csc_matrix(
        (np.array([1, 1, 1]), (np.array([0, 1, 2]), np.array([0, 0, 1])))
    )

    def sum_in_groups_side_effect(*args):
        return args[0] * Ih_table._csc_h_index_matrix

    Ih_table.sum_in_groups = sum_in_groups_side_effect
    return Ih_table


def mock_Ih_table():
    """A mock Ih table for testing a target function."""
    Ih_table = MagicMock()
    Ih_table.blocked_data_list = [mock_single_Ih_table()]
    return Ih_table


@pytest.fixture
def mock_restrained_component():
    """Mock a component with restraints."""
    component = Mock()
    component.n_params = 3
    component.calculate_restraints.return_value = [
        flex.double([1.0, 2.0, 3.0]),
        flex.double([0.1, 0.2, 0.3]),
    ]
    jacobian_restr = sparse.matrix(component.n_params, component.n_params)
    jacobian_restr[0, 0] = 1.0
    component.calculate_jacobian_restraints.return_value = [
        flex.double([1.0, 1.1, 1.2]),
        jacobian_restr,
        flex.double([1.0, 2.0, 3.0]),
    ]
    return component


@pytest.fixture
def mock_unrestrained_component():
    """Mock a component without restraints."""
    component = Mock()
    component.n_params = 5
    component.calculate_restraints.return_value = None
    component.calculate_jacobian_restraints.return_value = None
    return component


def _component_to_apm(component):
    mock_single_apm = Mock()
    mock_single_apm.components = {
        "restrained": {
            "object": component,
            "n_params": component.n_params,
            "start_idx": 0,
        }
    }
    mock_single_apm.n_active_params = component.n_params
    return mock_single_apm


def _single_to_multi_apm(single_apm):
    apm = MagicMock()
    apm.apm_list = [single_apm]
    apm.n_active_params = single_apm.n_active_params
    apm.apm_data = {0: {"start_idx": 0, "end_idx": single_apm.n_active_params}}
    return apm


@pytest.fixture
def mock_apm_restrained(mock_restrained_component):
    """Make a parameter manager with a restrained component"""
    return _single_to_multi_apm(_component_to_apm(mock_restrained_component))


@pytest.fixture
def mock_apm_unrestrained(mock_unrestrained_component):
    """Make a parameter manager with an unrestrained component"""
    return _single_to_multi_apm(_component_to_apm(mock_unrestrained_component))


def test_target_function_methods():
    """Test for the target methods required for the refinement engine."""
    target = ScalingTarget()
    r, w = target.compute_residuals(mock_single_Ih_table())
    assert r.size() == w.size()
    assert r == pytest.approx([-1.0, 0.0, 1.0])
    assert w == pytest.approx([1.0, 1.0, 1.0])

    f, g = target.compute_functional_gradients(mock_single_Ih_table())
    assert f == pytest.approx(2.0)
    assert g == pytest.approx([-19.04072398])

    r2, j, w2 = target.compute_residuals_and_gradients(mock_single_Ih_table())
    assert r == r2
    assert w == w2
    assert j.n_cols == 1 and j.n_rows == 3


def test_target_function_restraints_methods(mock_apm_restrained, mock_apm_unrestrained):
    """Test for the target restraints methods required for the refinement engine."""
    target = ScalingTarget()

    restraints = target.compute_restraints_residuals_and_gradients(mock_apm_restrained)
    assert len(restraints) == 3
    assert target.param_restraints is True
    assert list(restraints[0]) == [1.0, 1.1, 1.2]
    assert list(restraints[2]) == [1.0, 2.0, 3.0]

    restraints = target.compute_restraints_functional_gradients(mock_apm_restrained)
    assert len(restraints) == 2
    assert restraints[0] == 6.0
    assert list(restraints[1]) == [0.1, 0.2, 0.3]

    achieved = target.achieved()
    assert isinstance(achieved, bool)

    restraints = target.compute_restraints_residuals_and_gradients(
        mock_apm_unrestrained
    )
    assert restraints is None
    assert target.param_restraints is False
    # Check that if called again, returns None without doing calculations
    restraints = target.compute_restraints_residuals_and_gradients([])
    assert restraints is None

    target = ScalingTarget()  # Need to make new instance or won't calc restr as
    # param_restraints is set to False
    assert target.param_restraints is True
    restraints = target.compute_restraints_functional_gradients(mock_apm_unrestrained)
    assert restraints is None
    assert target.param_restraints is False


def test_target_rmsd_calculation(mock_apm_restrained, mock_apm_unrestrained):
    """Test the RMSD calculation, with and without restraints."""
    target = ScalingTarget()
    assert target.param_restraints is True
    # with input, expect residuals of [-1, 0, 1], weights of [1, 1, 1],
    # restraints of [1, 2, 3], so expect residual of sqrt((2+6)/3)
    rmsds = target.rmsds(mock_Ih_table(), mock_apm_restrained)
    assert len(rmsds) == 1
    assert rmsds[0] == pytest.approx((8.0 / 3.0) ** 0.5, abs=1e-6)
    assert target.param_restraints is True

    # test rmsd calculation without restraints
    rmsds = target.rmsds(mock_Ih_table(), mock_apm_unrestrained)
    assert len(rmsds) == 1
    assert rmsds[0] == pytest.approx((2.0 / 3.0) ** 0.5, abs=1e-6)
    assert target.param_restraints is False


def test_target_fixedIh():
    """Test the target function for targeted scaling (where Ih is fixed)."""

    target = ScalingTargetFixedIH()
    Ih_table = mock_Ih_table().blocked_data_list[0]
    R, _ = target.compute_residuals(Ih_table)
    expected_residuals = flex.double([-1.0, 0.0, 1.0])
    assert list(R) == pytest.approx(list(expected_residuals))
    _, G = target.compute_functional_gradients(Ih_table)
    assert list(G) == pytest.approx([-44.0])
    # Add in finite difference check

    Ih_table = mock_Ih_table().blocked_data_list[0]
    J = target.calculate_jacobian(Ih_table)
    assert J.n_cols == 1
    assert J.n_rows == 3
    assert J.non_zeroes == 3
    assert J[0, 0] == pytest.approx(-11.0)
    assert J[1, 0] == pytest.approx(-22.0)
    assert J[2, 0] == pytest.approx(-33.0)

    expected_rmsd = (
        flex.sum(flex.pow2(expected_residuals)) / len(expected_residuals)
    ) ** 0.5
    assert target._rmsds is None
    target.param_restraints = False  # don't try to use apm to get restraints
    assert target.rmsds(mock_Ih_table(), [])
    assert target._rmsds == pytest.approx([expected_rmsd])


# For testing the targetfunction calculations using finite difference methods,
# need to initialise real instances of the scaling datastructures to allow
# variation of the parameters and updating of the linked datastructures.


def test_target_gradient_calculation_finite_difference(
    small_reflection_table, single_exp, physical_param
):
    """Test the calculated gradients against a finite difference calculation."""
    model = PhysicalScalingModel.from_data(
        physical_param, single_exp, small_reflection_table
    )

    # need to 'add_data'
    model.configure_components(small_reflection_table, single_exp, physical_param)
    model.components["scale"].update_reflection_data()
    model.components["decay"].update_reflection_data()
    apm = multi_active_parameter_manager(
        ScalingTarget(),
        [model.components],
        [["scale", "decay"]],
        scaling_active_parameter_manager,
    )
    model.components["scale"].inverse_scales = flex.double([2.0, 1.0, 2.0])
    model.components["decay"].inverse_scales = flex.double([1.0, 1.0, 0.4])

    Ih_table = IhTable([small_reflection_table], single_exp.crystal.get_space_group())

    with patch.object(SingleScaler, "__init__", lambda x, y, z, k: None):
        scaler = SingleScaler(None, None, None)
        scaler._Ih_table = Ih_table

        # Now do finite difference check.
        target = ScalingTarget()

        scaler.update_for_minimisation(apm, 0)
        grad = target.calculate_gradients(scaler.Ih_table.blocked_data_list[0])
        res, _ = target.compute_residuals(scaler.Ih_table.blocked_data_list[0])

        assert (
            res > 1e-8
        ), """residual should not be zero, or the gradient test
        below will not really be working!"""

        # Now compare to finite difference
        f_d_grad = calculate_gradient_fd(target, scaler, apm)
        print(list(f_d_grad))
        print(list(grad))
        assert list(grad) == pytest.approx(list(f_d_grad))

        sel = f_d_grad > 1e-8
        assert sel, """assert sel has some elements, as finite difference grad should
        not all be zero, or the test will not really be working!
        (expect one to be zero for KB scaling example?)"""


def test_target_jacobian_calculation_finite_difference(
    physical_param, single_exp, large_reflection_table
):
    """Test the calculated jacobian against a finite difference calculation."""
    physical_param.physical.decay_correction = False
    model = PhysicalScalingModel.from_data(
        physical_param, single_exp, large_reflection_table
    )
    # need to 'add_data'
    model.configure_components(large_reflection_table, single_exp, physical_param)
    model.components["scale"].update_reflection_data()
    apm = multi_active_parameter_manager(
        ScalingTarget(),
        [model.components],
        [["scale"]],
        scaling_active_parameter_manager,
    )
    Ih_table = IhTable([large_reflection_table], single_exp.crystal.get_space_group())

    with patch.object(SingleScaler, "__init__", lambda x, y, z, k: None):
        scaler = SingleScaler(None, None, None)
        scaler._Ih_table = Ih_table

        target = ScalingTarget()
        scaler.update_for_minimisation(apm, 0)

        fd_jacobian = calculate_jacobian_fd(target, scaler, apm)
        r, jacobian, w = target.compute_residuals_and_gradients(
            scaler.Ih_table.blocked_data_list[0]
        )
        assert r == pytest.approx(
            [-50.0 / 3.0, 70.0 / 3.0, -20.0 / 3.0, 12.5, -2.5]
            + [-25.0, 0.0, -75.0, 0.0, 200.0]
        )
        assert w == pytest.approx(
            [0.1, 0.1, 0.1, 0.02, 0.1, 0.02, 0.01, 0.02, 0.01, 0.01]
        )

        n_rows = jacobian.n_rows
        n_cols = jacobian.n_cols

        print(jacobian)
        print(fd_jacobian)

        for i in range(0, n_rows):
            for j in range(0, n_cols):
                assert jacobian[i, j] == pytest.approx(fd_jacobian[i, j], abs=1e-4)


def calculate_gradient_fd(target, scaler, apm):
    """Calculate gradient array with finite difference approach."""
    delta = 1.0e-6
    gradients = flex.double([0.0] * apm.n_active_params)
    Ih_table = scaler.Ih_table.blocked_data_list[0]
    # iterate over parameters, varying one at a time and calculating the gradient
    for i in range(apm.n_active_params):
        new_x = copy.copy(apm.x)
        new_x[i] -= 0.5 * delta
        apm.set_param_vals(new_x)
        scaler.update_for_minimisation(apm, 0)
        R_low = np.square(target.calculate_residuals(Ih_table)) * Ih_table.weights
        new_x[i] += delta
        apm.set_param_vals(new_x)
        scaler.update_for_minimisation(apm, 0)
        R_upper = np.square(target.calculate_residuals(Ih_table)) * Ih_table.weights
        new_x[i] -= 0.5 * delta
        apm.set_param_vals(new_x)
        scaler.update_for_minimisation(apm, 0)
        gradients[i] = (np.sum(R_upper) - np.sum(R_low)) / delta
    return gradients


def calculate_jacobian_fd(target, scaler, apm, block_id=0):
    """Calculate jacobian matrix with finite difference approach."""
    delta = 1.0e-7
    jacobian = sparse.matrix(
        scaler.Ih_table.blocked_data_list[block_id].size, apm.n_active_params
    )
    Ih_table = scaler.Ih_table.blocked_data_list[block_id]
    # iterate over parameters, varying one at a time and calculating the residuals
    for i in range(apm.n_active_params):
        new_x = copy.copy(apm.x)
        new_x[i] -= 0.5 * delta
        apm.set_param_vals(new_x)
        scaler.update_for_minimisation(apm, 0)
        R_low = target.calculate_residuals(Ih_table)  # unweighted unsquared residual
        new_x[i] += delta
        apm.set_param_vals(new_x)
        scaler.update_for_minimisation(apm, 0)
        R_upper = target.calculate_residuals(Ih_table)  # unweighted unsquared residual
        new_x[i] -= 0.5 * delta
        apm.set_param_vals(new_x)
        scaler.update_for_minimisation(apm, 0)
        fin_difference = (R_upper - R_low) / delta
        for j in range(fin_difference.size):
            jacobian[j, i] = fin_difference[j]
    return jacobian
