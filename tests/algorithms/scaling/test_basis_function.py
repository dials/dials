"""
Test for the basis function module.
"""

from __future__ import annotations

import pytest

from dials.algorithms.scaling.basis_functions import RefinerCalculator
from dials.algorithms.scaling.model.components.scale_components import (
    SingleBScaleFactor,
    SingleScaleFactor,
)
from dials.algorithms.scaling.parameter_handler import scaling_active_parameter_manager
from dials.array_family import flex


@pytest.fixture
def small_reflection_table():
    """Generate reflection table to test the basis function."""
    reflections = flex.reflection_table()
    reflections["d"] = flex.double([2.0, 0.8, 2.0])  # don't change
    reflections["xyzobs.px.value"] = flex.vec3_double(
        [(0.0, 0.0, 0.0), (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)]
    )
    reflections["id"] = flex.int(3, 0)
    return reflections


def test_RefinerCalculator(small_reflection_table):
    """Test for the RefinerCalculator class. This calculates scale factors and
    derivatives for reflections based on the model components."""

    # To test the basis function, need a scaling active parameter manager - to set
    # this up we need a components dictionary with some reflection data.

    # Let's use KB model components for simplicity - and have an extra fake 'abs'
    # component.
    rt = small_reflection_table
    components = {
        "scale": SingleScaleFactor(flex.double([1.0])),
        "decay": SingleBScaleFactor(flex.double([0.0])),
        "abs": SingleScaleFactor(flex.double([1.0])),
    }  # Create empty components.
    components["scale"].data = {"id": rt["id"]}
    components["decay"].data = {"d": rt["d"]}
    components["abs"].data = {"id": rt["id"]}
    for component in components.values():
        component.update_reflection_data()  # Add some data to components.

    apm = scaling_active_parameter_manager(components, ["decay", "scale"])

    # First test that scale factors can be successfully updated.
    # Manually change the parameters in the apm.
    decay = components["decay"]  # Define alias
    _ = components["scale"]  # Define alias
    # Note, order of params in apm.x depends on order in scaling model components.
    new_B = 1.0
    new_S = 2.0
    apm.set_param_vals(flex.double([new_S, new_B]))
    s, d = RefinerCalculator.calculate_scales_and_derivatives(apm, 0)
    slist, dlist = RefinerCalculator._calc_component_scales_derivatives(apm, 0)
    # Now test that the inverse scale factor is correctly calculated.
    calculated_sfs = s
    assert list(calculated_sfs) == pytest.approx(
        list(new_S * flex.exp(new_B / (2.0 * flex.pow2(decay.d_values[0]))))
    )

    # Now check that the derivative matrix is correctly calculated.
    calc_derivs = d
    assert calc_derivs[0, 0] == dlist[0][0, 0] * slist[1][0]
    assert calc_derivs[1, 0] == dlist[0][1, 0] * slist[1][1]
    assert calc_derivs[2, 0] == dlist[0][2, 0] * slist[1][2]
    assert calc_derivs[0, 1] == dlist[1][0, 0] * slist[0][0]
    assert calc_derivs[1, 1] == dlist[1][1, 0] * slist[0][1]
    assert calc_derivs[2, 1] == dlist[1][2, 0] * slist[0][2]

    # Repeat the test when there is only one active parameter.
    # First reset the parameters
    components["decay"].parameters = flex.double([0.0])
    components["scale"].parameters = flex.double([1.0])
    components["abs"].parameters = flex.double([1.0])
    components["decay"].calculate_scales_and_derivatives()
    components["scale"].calculate_scales_and_derivatives()
    components["abs"].calculate_scales_and_derivatives()

    # Now generate a parameter manager for a single component.
    apm = scaling_active_parameter_manager(components, ["scale"])
    new_S = 2.0
    apm.set_param_vals(flex.double(components["scale"].n_params, new_S))
    s, d = RefinerCalculator.calculate_scales_and_derivatives(apm, 0)
    slist, dlist = RefinerCalculator._calc_component_scales_derivatives(apm, 0)
    # Test that the scales and derivatives were correctly calculated
    assert list(s) == list([new_S] * slist[0].size())
    assert d[0, 0] == dlist[0][0, 0]
    assert d[1, 0] == dlist[0][1, 0]
    assert d[2, 0] == dlist[0][2, 0]

    # Test again for two components.
    components["decay"].parameters = flex.double([0.0])
    components["scale"].parameters = flex.double([1.0])
    components["abs"].parameters = flex.double([1.0])
    components["decay"].calculate_scales_and_derivatives()
    components["scale"].calculate_scales_and_derivatives()
    components["abs"].calculate_scales_and_derivatives()

    apm = scaling_active_parameter_manager(components, ["scale", "decay"])
    _, __ = RefinerCalculator.calculate_scales_and_derivatives(apm, 0)

    # Test for no components
    apm = scaling_active_parameter_manager(components, [])
    _, d = RefinerCalculator.calculate_scales_and_derivatives(apm, 0)
    assert d.n_cols == 0 and d.n_rows == 0
