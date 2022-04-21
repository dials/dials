from __future__ import annotations

from dials.array_family.flex import reflection_table
from dials.model.data import MultiPanelImageVolume
from dials_algorithms_integration_sum_ext import *  # noqa: F403; lgtm
from dials_algorithms_integration_sum_ext import sum_image_volume

__all__ = (  # noqa: F405
    "SummationDouble",
    "SummationFloat",
    "integrate_by_summation",
    "sum_image_volume",
    "sum_integrate_and_update_table",
)


def sum_integrate_and_update_table(
    reflections: reflection_table, image_volume: MultiPanelImageVolume = None
) -> reflection_table:
    """Perform 3D summation integration and update a reflection table.

    Arguments:
        reflections: The reflections to integrate

    Returns:
        The integrated reflections
    """

    # Integrate and return the reflections
    if image_volume is None:
        intensity = reflections["shoebox"].summed_intensity()
    else:
        intensity = sum_image_volume(reflections, image_volume)
    reflections["intensity.sum.value"] = intensity.observed_value()
    reflections["intensity.sum.variance"] = intensity.observed_variance()
    reflections["background.sum.value"] = intensity.background_value()
    reflections["background.sum.variance"] = intensity.background_variance()
    success = intensity.observed_success()
    reflections.set_flags(success, reflections.flags.integrated_sum)
    return success
