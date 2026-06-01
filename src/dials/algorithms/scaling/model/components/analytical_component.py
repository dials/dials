from __future__ import annotations

from scitbx import sparse
from scitbx.array_family import flex

from dials.algorithms.scaling.model.components.scale_components import (
    ScaleComponentBase,
)


class AnalyticalComponent(ScaleComponentBase):
    null_parameter_value = 1.0

    def __init__(self, initial_values, parameter_esds=None):
        """Set the initial parameter values, parameter esds and n_params."""
        super().__init__(initial_values, parameter_esds)
        self._analytical_corrections = []

    @ScaleComponentBase.data.setter
    def data(self, data):
        """Set the data dict in the parent class."""
        assert set(data.keys()) == {"correction"}, set(data.keys())
        self._data = data

    @property
    def analytical_corrections(self):
        return self._analytical_corrections

    @property
    def free_parameters(self):
        return flex.double([])

    @free_parameters.setter
    def free_parameters(self, parameters):
        self._parameters = flex.double([])

    @property
    def free_parameter_esds(self):
        """Return the estimated standard deviations of the parameters."""
        return flex.double([])

    @free_parameter_esds.setter
    def free_parameter_esds(self, esds):
        self._parameter_esds = flex.double([])

    def update_reflection_data(self, selection=None, block_selections=None):
        data = self.data["correction"]
        if selection:
            self._analytical_corrections = [data.select(selection)]
        elif block_selections:
            self._analytical_corrections = [
                data.select(sel) for sel in block_selections
            ]
        else:
            self._analytical_corrections = [data]

        self._n_refl = [s.size() for s in self._analytical_corrections]

    def calculate_scales_and_derivatives(self, block_id=0):
        """Calculate and return inverse scales and derivatives for a given block."""
        derivatives = sparse.matrix(self.n_refl[block_id], 0)
        return self._analytical_corrections[block_id], derivatives

    def calculate_scales(self, block_id=0):
        """Calculate and return inverse scales for a given block."""
        return self._analytical_corrections[block_id]
