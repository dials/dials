from __future__ import annotations

from cctbx import uctbx
from scitbx.array_family import flex


def quasi_normalisation(intensities):
    """Quasi-normalisation of the input intensities.

    Args:
        intensities (cctbx.miller.array): The intensities to be normalised.

    Returns:
        cctbx.miller.array: The normalised intensities.
    """
    # handle negative reflections to minimise effect on mean I values.
    work = intensities.deep_copy()
    work.data().set_selected(work.data() < 0.0, 0.0)

    # set up binning objects
    if work.size() > 20000:
        n_refl_shells = 20
    elif work.size() > 15000:
        n_refl_shells = 15
    else:
        n_refl_shells = 10
    d_star_sq = work.d_star_sq().data()
    d_star_sq_max = flex.max(d_star_sq)
    d_star_sq_min = flex.min(d_star_sq)
    span = d_star_sq_max - d_star_sq_min
    d_star_sq_max += span * 1e-6
    d_star_sq_min -= span * 1e-6
    d_star_sq_step = (d_star_sq_max - d_star_sq_min) / n_refl_shells
    work.setup_binner_d_star_sq_step(
        d_min=uctbx.d_star_sq_as_d(d_star_sq_min),  # cctbx/cctbx_project#588
        d_max=uctbx.d_star_sq_as_d(d_star_sq_max),  # cctbx/cctbx_project#588
        d_star_sq_step=d_star_sq_step,
        auto_binning=False,
    )
    normalisations = work.intensity_quasi_normalisations()
    return intensities.customized_copy(
        data=(intensities.data() / normalisations.data()),
        sigmas=(intensities.sigmas() / normalisations.data()),
    )
