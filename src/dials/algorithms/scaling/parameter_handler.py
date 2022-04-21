"""
Extension to general active parameter manager for scaling and function
to use a scaler to determine the correct call to the apm factories.
"""

from __future__ import annotations

from dials.algorithms.scaling.active_parameter_managers import (
    ParameterManagerGenerator,
    active_parameter_manager,
)


class scaling_active_parameter_manager(active_parameter_manager):
    """
    Adds an extra property to the apm to avoid a repetitive calculation during
    mimimisation cycles for scaling.
    """

    def __init__(self, components, selection_list):
        self.constant_g_values = None
        for component, obj in components.items():
            if component not in selection_list:
                n_blocks = len(obj.n_refl)
                if self.constant_g_values is None:
                    self.constant_g_values = [None] * n_blocks
                    for n in range(n_blocks):
                        self.constant_g_values[n] = obj.calculate_scales(n)
                else:
                    for n in range(n_blocks):
                        self.constant_g_values[n] *= obj.calculate_scales(n)
        super().__init__(components, selection_list)
        n_obs = []
        for component in components:
            obs_in_component = []
            for n_refl in components[component].n_refl:
                obs_in_component.append(n_refl)
            n_obs.append(obs_in_component)
        assert all(i == n_obs[0] for i in n_obs)
        n_obs = []
        for component in components:
            n_obs.append(components[component].n_refl)
        self.n_obs = n_obs[0]  # list of length n_blocks


class ScalingParameterManagerGenerator(ParameterManagerGenerator):

    """Class to generate parameter manager for scaling."""

    def __init__(self, data_managers, target, mode, shared=None):
        super().__init__(
            data_managers, scaling_active_parameter_manager, target, mode, shared
        )
