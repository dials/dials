"""Scaling model components for time-of-flight data."""

from __future__ import annotations

from collections.abc import Iterable

from scitbx import sparse

from dials.algorithms.scaling.model.components.scale_components import (
    ScaleComponentBase,
)
from dials.array_family import flex


class PanelScaleComponent(ScaleComponentBase):
    """One scale factor per detector panel.

    The inverse scale factor of a reflection is the parameter of the panel it
    was recorded on. Reference panel and those without reflections held fixed.
    """

    null_parameter_value = 1.0

    def __init__(
        self,
        initial_values: flex.double,
        parameter_esds: flex.double | None = None,
        fixed_panels: Iterable[int] = (),
    ) -> None:
        super().__init__(initial_values, parameter_esds)
        self._fixed_panels = sorted(set(fixed_panels))
        assert all(0 <= p < self.n_params for p in self._fixed_panels)
        self._free = flex.bool(self.n_params, True)
        for p in self._fixed_panels:
            self._free[p] = False
        # column of each panel in the derivatives matrix (-1 if fixed)
        self._column = []
        n = 0
        for free in self._free:
            self._column.append(n if free else -1)
            n += int(free)
        self._panels = []

    @property
    def fixed_panels(self) -> list[int]:
        return self._fixed_panels

    @ScaleComponentBase.data.setter
    def data(self, data: dict[str, flex.size_t]) -> None:
        assert set(data.keys()) == {"panel"}, set(data.keys())
        self._data = data

    @property
    def free_parameters(self) -> flex.double:
        return self._parameters.select(self._free)

    @free_parameters.setter
    def free_parameters(self, parameters: flex.double) -> None:
        self._parameters.set_selected(self._free, parameters)

    @property
    def free_parameter_esds(self) -> flex.double | None:
        if self._parameter_esds is None:
            return None
        return self._parameter_esds.select(self._free)

    @free_parameter_esds.setter
    def free_parameter_esds(self, esds: flex.double) -> None:
        assert len(esds) == len(self.free_parameters)
        if self._parameter_esds is None:
            self._parameter_esds = flex.double(self.n_params, 0.0)
        self._parameter_esds.set_selected(self._free, esds)

    def update_reflection_data(
        self,
        selection: flex.bool | flex.size_t | None = None,
        block_selections: list[flex.size_t] | None = None,
    ) -> None:
        panels = flex.size_t(list(self.data["panel"]))
        if selection:
            self._panels = [panels.select(selection)]
        elif block_selections:
            self._panels = [panels.select(sel) for sel in block_selections]
        else:
            self._panels = [panels]
        self._n_refl = [p.size() for p in self._panels]

    def calculate_scales(self, block_id: int = 0) -> flex.double:
        return self._parameters.select(self._panels[block_id])

    def calculate_scales_and_derivatives(
        self, block_id: int = 0
    ) -> tuple[flex.double, sparse.matrix]:
        panels = self._panels[block_id]
        derivatives = sparse.matrix(panels.size(), len(self.free_parameters))
        for i, p in enumerate(panels):
            col = self._column[p]
            if col >= 0:
                derivatives[i, col] = 1.0
        return self.calculate_scales(block_id), derivatives
