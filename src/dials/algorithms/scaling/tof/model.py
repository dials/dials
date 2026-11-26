"""
Scaling model for time-of-flight Laue data.

Wavelength dependence of the intensities is expected to have been removed by
dials.lambda_norm, so datasets are treated as effectively monochromatic.

- scale: one scale factor per run
- decay: one relative B-factor per run (off by default)
- absorption: a spherical-harmonic absorption surface in the crystal frame,
  shared across all runs and scaled linearly with wavelength (off by default).
  This accounts for the absorption not covered the simple spherical geometrical
  absorption correction in dials.tof_itegrate (e.g. non-spherical crystals)
- panel: one scale factor per detector panel, shared across all runs (off by
  default). One panel is a fixed reference, as an overall scale is already
  given by the run scales.
- analytical: a fixed correction from the "analytical_correction" column
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from libtbx import phil
from scitbx import sparse

from dials.algorithms.scaling.model.components.analytical_component import (
    AnalyticalComponent,
)
from dials.algorithms.scaling.model.components.scale_components import (
    SHScaleComponent,
    SingleBScaleFactor,
    SingleScaleFactor,
)
from dials.algorithms.scaling.model.model import (
    ScalingModelBase,
    base_model_phil_str,
)
from dials.algorithms.scaling.plots import (
    plot_absorption_parameters,
    plot_absorption_plots,
)
from dials.algorithms.scaling.scaling_utilities import sph_harm_table
from dials.algorithms.scaling.tof.components import PanelScaleComponent
from dials.array_family import flex

if TYPE_CHECKING:
    from dxtbx.model import Experiment, ExperimentList

logger = logging.getLogger("dials")

WAVELENGTH = "wavelength_cal"

tof_model_phil_str = (
    """\
decay_correction = False
    .type = bool
    .help = "Refine one relative B-factor per run."
    .expert_level = 1
absorption_correction = False
    .type = bool
    .help = "Refine a spherical-harmonic absorption surface in the crystal"
            "frame, shared across all runs. Requires more than one run."
    .expert_level = 1
absorption_wavelength_dependence = *linear none
    .type = choice
    .help = "linear: the absorption surface is scaled by lambda/lambda_ref,"
            "following the 1/v law for neutron absorption cross-sections."
            "none: the same surface is applied at every wavelength."
    .expert_level = 1
lmax = 4
    .type = int(value_min=2)
    .help = "Number of spherical harmonics to include for absorption correction."
    .expert_level = 2
surface_weight = 1e3
    .type = float(value_min=0.0)
    .help = "Restraint weight applied to spherical harmonic terms in the"
            "absorption correction."
    .expert_level = 2
reference_wavelength = None
    .type = float(value_min=0)
    .help = "Wavelength (A) at which the absorption surface is defined."
            "Default: the median wavelength of all observations."
    .expert_level = 2
panel_correction = False
    .type = bool
    .help = "Refine one scale factor per detector panel, shared across all"
            "runs. Requires more than one run."
    .expert_level = 1
panel_reference = None
    .type = int(value_min=0)
    .help = "The panel whose scale is fixed at 1. Default: the panel with the"
            "most observations."
    .expert_level = 2
panel_fixed = None
    .type = ints(value_min=0)
    .help = "Panels whose scale is held at 1, as well as the reference panel."
            "Panels without observations are added automatically."
    .expert_level = 2
require_lambda_norm = True
    .type = bool
    .help = "Require the data to have been wavelength-normalised with"
            "dials.lambda_norm."
    .expert_level = 1
"""
    + base_model_phil_str
)


def _wavelength_scaled_sph_harm_table(
    reflection_table: flex.reflection_table, lmax: int, lambda_ref: float
) -> sparse.matrix:
    """The spherical harmonic table with each reflection's column scaled by
    lambda/lambda_ref, so that the standard SH component evaluates
    S = 1 + (lambda/lambda_ref) * sum_lm Clm Ylm."""
    table = sph_harm_table(reflection_table, lmax)
    factors = reflection_table[WAVELENGTH] / lambda_ref
    n = factors.size()
    diagonal = sparse.matrix(n, n)
    for i, f in enumerate(factors):
        diagonal[i, i] = f
    return table * diagonal


class TOFScalingModel(ScalingModelBase):
    """Scaling model for time-of-flight Laue data"""

    id_ = "tof"

    phil_scope = phil.parse(tof_model_phil_str)

    def __init__(
        self,
        parameters_dict: dict[str, dict[str, Any]],
        configdict: dict[str, Any],
        is_scaled: bool = False,
    ) -> None:
        super().__init__(configdict, is_scaled)
        if "scale" in configdict["corrections"]:
            self._components["scale"] = SingleScaleFactor(
                parameters_dict["scale"]["parameters"],
                parameters_dict["scale"]["parameter_esds"],
            )
        if "decay" in configdict["corrections"]:
            self._components["decay"] = SingleBScaleFactor(
                parameters_dict["decay"]["parameters"],
                parameters_dict["decay"]["parameter_esds"],
            )
        if "absorption" in configdict["corrections"]:
            self._components["absorption"] = SHScaleComponent(
                parameters_dict["absorption"]["parameters"],
                parameters_dict["absorption"]["parameter_esds"],
            )
        if "panel" in configdict["corrections"]:
            self._components["panel"] = PanelScaleComponent(
                parameters_dict["panel"]["parameters"],
                parameters_dict["panel"]["parameter_esds"],
                fixed_panels=configdict["fixed_panels"],
            )
        if "analytical" in configdict["corrections"]:
            self._components["analytical"] = AnalyticalComponent(flex.double([]))

    @property
    def consecutive_refinement_order(self) -> list[list[str]]:
        """:obj:`list`: a nested list of component names to indicate scaling order."""
        return [["scale", "decay", "panel"], ["absorption"]]

    def configure_components(
        self,
        reflection_table: flex.reflection_table,
        experiment: Experiment,
        params: phil.scope_extract,
    ) -> None:
        """Add the required reflection table data to the model components."""
        if "scale" in self.components:
            self.components["scale"].data = {"id": reflection_table["id"]}
        if "decay" in self.components:
            self.components["decay"].data = {"d": reflection_table["d"]}
        if "absorption" in self.components:
            self._configure_absorption(reflection_table)
        if "panel" in self.components:
            self.components["panel"].data = {"panel": reflection_table["panel"]}
        if "analytical" in self.components:
            if "analytical_correction" not in reflection_table:
                raise ValueError(
                    "The reflection table must contain the column 'analytical_correction'"
                    + "\nwhen using the option analytical_correction=True."
                )
            self.components["analytical"].data = {
                "correction": reflection_table["analytical_correction"]
            }

    def _configure_absorption(self, reflection_table: flex.reflection_table) -> None:
        if "s0c" not in reflection_table or "s1c" not in reflection_table:
            raise ValueError(
                "Crystal-frame vectors (s0c, s1c) are needed for the ToF absorption"
                " correction"
            )
        lmax = self._configdict["lmax"]
        if self._configdict["absorption_wavelength_dependence"] == "linear":
            table = _wavelength_scaled_sph_harm_table(
                reflection_table, lmax, self._configdict["reference_wavelength"]
            )
        else:
            table = sph_harm_table(reflection_table, lmax)
        component = self.components["absorption"]
        component.data = {"sph_harm_table": table}
        restraints = flex.double([])
        for l in range(1, lmax + 1):
            restraints.extend(flex.double([1.0] * ((2 * l) + 1)))
        restraints *= self._configdict["abs_surface_weight"]
        component.parameter_restraints = restraints

    def get_shared_components(self) -> tuple[str, ...] | None:
        shared = tuple(
            name
            for name in ("absorption", "panel")
            if name in self.components and name not in self.fixed_components
        )
        return shared or None

    def update(self, params: phil.scope_extract) -> None:
        """Update the model if new options chosen in the phil scope."""
        if params.tof.correction.fix:
            self.fixed_components = params.tof.correction.fix
        if "absorption" in self.components:
            self._configdict["abs_surface_weight"] = params.tof.surface_weight

    def plot_model_components(
        self, reflection_table: flex.reflection_table | None = None
    ) -> dict[str, dict]:
        d = {}
        if "panel" in self.components:
            d.update(plot_panel_scales(self))
        if "absorption" in self.components:
            d.update(plot_absorption_parameters(self))
            d.update(plot_absorption_plots(self, reflection_table))
        return d

    @classmethod
    def from_data(
        cls,
        params: phil.scope_extract,
        experiment: Experiment,
        reflection_table: flex.reflection_table,
    ) -> TOFScalingModel:
        params = params.tof
        configdict = {"corrections": ["scale"]}
        parameters_dict = {
            "scale": {"parameters": flex.double([1.0]), "parameter_esds": None}
        }
        if params.decay_correction:
            configdict["corrections"].append("decay")
            parameters_dict["decay"] = {
                "parameters": flex.double([0.0]),
                "parameter_esds": None,
            }
        if params.absorption_correction:
            configdict["corrections"].append("absorption")
            lmax = params.lmax
            wavelength_dependence = params.absorption_wavelength_dependence
            configdict.update(
                {
                    "lmax": lmax,
                    "abs_surface_weight": params.surface_weight,
                    "absorption_wavelength_dependence": wavelength_dependence,
                }
            )
            if wavelength_dependence == "linear":
                if not params.reference_wavelength:
                    raise ValueError(
                        "A reference wavelength is needed for the wavelength-dependent"
                        " absorption correction"
                    )
                configdict["reference_wavelength"] = params.reference_wavelength
            n_abs_param = (2 * lmax) + (lmax**2)
            parameters_dict["absorption"] = {
                "parameters": flex.double(n_abs_param, 0.0),
                "parameter_esds": None,
            }
        if params.panel_correction:
            n_panels = len(experiment.detector)
            if params.panel_reference is None:
                raise ValueError("A reference panel is needed for the panel correction")
            if params.panel_reference >= n_panels:
                raise ValueError(
                    f"panel_reference={params.panel_reference}, but the detector"
                    f" has {n_panels} panels"
                )
            configdict["corrections"].append("panel")
            configdict["fixed_panels"] = sorted(
                {params.panel_reference, *(params.panel_fixed or [])}
            )
            parameters_dict["panel"] = {
                "parameters": flex.double(n_panels, 1.0),
                "parameter_esds": None,
            }
        if params.analytical_correction:
            if "analytical_correction" not in reflection_table:
                raise ValueError(
                    "The reflection table must contain the column 'analytical_correction'"
                    + "\nwhen using the option analytical_correction=True."
                )
            configdict["corrections"].append("analytical")

        model = cls(parameters_dict, configdict)
        if params.correction.fix:
            model.fixed_components = params.correction.fix
        return model

    @classmethod
    def from_dict(cls, obj: dict[str, Any]) -> TOFScalingModel:
        if obj["__id__"] != cls.id_:
            raise RuntimeError(f"expected __id__ {cls.id_}, got {obj['__id__']}")
        configdict = obj["configuration_parameters"]
        parameters_dict = {}
        for name in ("scale", "decay", "absorption", "panel"):
            if name in configdict["corrections"]:
                esds = None
                if "est_standard_devs" in obj[name]:
                    esds = flex.double(obj[name]["est_standard_devs"])
                parameters_dict[name] = {
                    "parameters": flex.double(obj[name]["parameters"]),
                    "parameter_esds": esds,
                }
        return cls(parameters_dict, configdict, is_scaled=True)


def plot_panel_scales(model: TOFScalingModel) -> dict[str, dict]:
    """A bar chart of the panel scale factors."""
    component = model.components["panel"]
    params = list(component.parameters)
    bar = {
        "x": list(range(len(params))),
        "y": [round(v, 5) for v in params],
        "type": "bar",
        "name": "panel scale",
    }
    if component.parameter_esds:
        bar["error_y"] = {
            "type": "data",
            "array": [round(v, 5) for v in component.parameter_esds],
        }
    return {
        "panel_scales": {
            "data": [bar],
            "layout": {
                "title": "Panel scale factors (shared across runs)",
                "xaxis": {"title": "Panel", "dtick": 1},
                "yaxis": {"title": "Inverse scale factor"},
            },
            "help": "One inverse scale factor per detector panel, shared by all"
            " runs. Fixed panels (the reference, and any panel without"
            f" observations) are held at 1: {component.fixed_panels}.",
        }
    }


def check_panel_scales(
    experiments: ExperimentList, max_ratio: float = 2.0
) -> list[int]:
    """Warn about panels whose refined scale is far from that of the others,
    which usually means bad data on that panel.
    Returns the list of bad panels."""
    models = [e.scaling_model for e in experiments if e.scaling_model]
    models = [m for m in models if m.id_ == "tof" and "panel" in m.components]
    if not models:
        return []
    scales = list(models[0].components["panel"].parameters)
    median = sorted(scales)[len(scales) // 2]
    outliers = [
        i for i, s in enumerate(scales) if not (1 / max_ratio < s / median < max_ratio)
    ]
    if outliers:
        logger.warning(
            "Warning: the scale of panel(s) %s differs from the median panel scale"
            " by more than a factor of %s: %s. This often means poor data on the"
            " panel (e.g. mis-predicted reflections), which the panel scale will"
            " amplify; consider excluding it.",
            ", ".join(str(i) for i in outliers),
            max_ratio,
            ", ".join(f"{scales[i]:.3f}" for i in outliers),
        )
    return outliers
