"""
Preparation of time-of-flight Laue data for dials.scale.

Methods check the input is suitable, choose the ToF scaling model,
and add the crystal-frame vectors needed for an absorption correction.

The data are expected to have been wavelength-normalised by dials.lambda_norm,
so that scaling can treat them as monochromatic.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from dxtbx import flumpy
from dxtbx.model import ExperimentType
from libtbx import Auto

from dials.algorithms.scaling.scaling_utilities import align_axis_along_z
from dials.algorithms.scaling.tof.model import WAVELENGTH
from dials.array_family import flex

if TYPE_CHECKING:
    from dxtbx.model import Experiment, ExperimentList
    from libtbx.phil import scope_extract

logger = logging.getLogger("dials")

LAMBDA_NORM_COLUMN = "lambda_norm.scale"
TOF_COMPATIBLE_MODELS = ("tof", "KB")


def is_tof_experiment(experiment: Experiment) -> bool:
    return experiment.get_type() == ExperimentType.TOF


def prepare_tof_input(
    params: scope_extract,
    experiments: ExperimentList,
    reflections: list[flex.reflection_table],
) -> tuple[scope_extract, ExperimentList, list[flex.reflection_table]]:
    """Check that ToF experiments can be scaled, and prepare them.

    Must be called after the reflections have been split into one table per
    experiment.
    """
    if not all(is_tof_experiment(e) for e in experiments):
        raise ValueError(
            "Time-of-flight and other experiment types cannot be scaled together"
        )
    logger.info("Time-of-flight data: preparing for ToF scaling.")
    _check_options(params, experiments)
    _check_reflections(params, reflections)
    if params.tof.absorption_correction:
        if params.tof.absorption_wavelength_dependence == "linear":
            _set_reference_wavelength(params, reflections)
    if params.tof.panel_correction:
        _set_panel_options(params, experiments, reflections)
    if params.tof.absorption_correction or _existing_models_have_absorption(
        params, experiments
    ):
        for expt, table in zip(experiments, reflections):
            add_tof_crystal_frame_vectors(table, expt)
    return params, experiments, reflections


def _existing_models_have_absorption(
    params: scope_extract, experiments: ExperimentList
) -> bool:
    """Whether models from a previous scaling run, that will be reused, have an
    absorption correction."""
    if params.overwrite_existing_models:
        return False
    return any(
        e.scaling_model and "absorption" in e.scaling_model.components
        for e in experiments
    )


def _check_options(params: scope_extract, experiments: ExperimentList) -> None:
    if params.model in (None, Auto, "auto", "Auto"):
        params.model = "tof"
        logger.info("Set model=tof for time-of-flight data")
    elif params.model not in TOF_COMPATIBLE_MODELS:
        raise ValueError(
            f"The {params.model} model cannot be used for time-of-flight data:"
            f" use one of {', '.join(TOF_COMPATIBLE_MODELS)}"
        )
    if len(experiments) == 1 and not params.scaling_options.reference:
        raise ValueError(
            "Time-of-flight scaling needs more than one run: with one run, the"
            " single scale factor has nothing to be refined against"
        )
    if params.tof.absorption_correction and len(experiments) == 1:
        raise ValueError(
            "The ToF absorption correction needs more than one run, as the beam"
            " direction is fixed in the crystal frame within a run"
        )
    if any(params.exclude_images):
        raise ValueError(
            "exclude_images cannot be used for time-of-flight data, where images"
            " are time-of-flight bins: use wavelength_range in dials.tof_integrate"
            " or dials.lambda_norm instead"
        )
    if params.output.merged_mtz:
        raise ValueError(
            "output.merged_mtz is not supported for time-of-flight data:"
            " use output.unmerged_mtz, or dials.export"
        )
    if (
        params.filtering.method == "deltacchalf"
        and params.filtering.deltacchalf.mode == "image_group"
    ):
        raise ValueError(
            "filtering.deltacchalf.mode=image_group cannot be used for"
            " time-of-flight data: use mode=dataset"
        )


def _check_reflections(
    params: scope_extract, reflections: list[flex.reflection_table]
) -> None:
    for table in reflections:
        if LAMBDA_NORM_COLUMN not in table:
            if params.tof.require_lambda_norm:
                raise ValueError(
                    "Time-of-flight data must be wavelength-normalised before scaling:"
                    " run dials.lambda_norm first (or set"
                    " tof.require_lambda_norm=False)"
                )
            logger.warning(
                "Warning: these time-of-flight data have not been wavelength-normalised"
            )
            break
    if (
        params.tof.absorption_correction
        and params.tof.absorption_wavelength_dependence == "linear"
    ):
        if not all(WAVELENGTH in table for table in reflections):
            raise ValueError(
                f"A {WAVELENGTH} column is needed for the wavelength-dependent"
                " absorption correction"
            )


def _set_reference_wavelength(
    params: scope_extract, reflections: list[flex.reflection_table]
) -> None:
    if params.tof.reference_wavelength:
        return
    wavelengths = flex.double()
    for table in reflections:
        wavelengths.extend(table[WAVELENGTH])
    params.tof.reference_wavelength = float(np.median(flumpy.to_numpy(wavelengths)))
    logger.info(
        "Absorption surface reference wavelength: %.3f Å",
        params.tof.reference_wavelength,
    )


def _set_panel_options(
    params: scope_extract,
    experiments: ExperimentList,
    reflections: list[flex.reflection_table],
) -> None:
    """Choose the reference panel, if not given, and fix panels without
    observations, so that every free panel scale is determined by the data."""
    n_panels = {len(e.detector) for e in experiments}
    if len(n_panels) != 1:
        raise ValueError(
            "The panel correction needs every run to have the same number of panels"
        )
    n_panels = n_panels.pop()
    if not all("panel" in table for table in reflections):
        raise ValueError("A panel column is needed for the panel correction")
    counts = np.zeros(n_panels, dtype=int)
    for table in reflections:
        good = ~table.get_flags(table.flags.bad_for_scaling, all=False)
        panels = flumpy.to_numpy(table["panel"].select(good)).astype(int)
        counts += np.bincount(panels, minlength=n_panels)[:n_panels]
    empty = [int(p) for p in np.flatnonzero(counts == 0)]
    if params.tof.panel_reference is None:
        params.tof.panel_reference = int(np.argmax(counts))
    elif counts[params.tof.panel_reference] == 0:
        raise ValueError(
            f"The reference panel {params.tof.panel_reference} has no observations"
        )
    params.tof.panel_fixed = sorted({*(params.tof.panel_fixed or []), *empty})
    logger.info(
        "Panel scale correction: reference panel %d (%d observations)%s",
        params.tof.panel_reference,
        counts[params.tof.panel_reference],
        f"; panels without observations, fixed: {empty}" if empty else "",
    )


def add_tof_crystal_frame_vectors(
    reflection_table: flex.reflection_table, experiment: Experiment
) -> flex.reflection_table:
    """Add the s0c and s1c columns: the unit beam and diffracted beam directions
    in the crystal frame.

    The crystal is stationary during a ToF run, and its orientation is given by
    the goniometer setting and fixed rotations, so there is no scan rotation.
    The frame is then rotated to put the lab x-axis along z, as for rotation data.
    """
    if experiment.goniometer:
        gonio = experiment.goniometer
        fixed_rotation = np.array(gonio.get_fixed_rotation()).reshape(3, 3)
        setting_rotation = np.array(gonio.get_setting_rotation()).reshape(3, 3)
        R = setting_rotation @ fixed_rotation
    else:
        R = np.identity(3)

    s0 = np.array(experiment.beam.get_sample_to_source_direction())
    s1 = flumpy.to_numpy(reflection_table["s1"])
    lengths = np.linalg.norm(s1, axis=1)
    non_zero = lengths > 0.0
    s0c = np.zeros((len(reflection_table), 3))
    s1c = np.zeros((len(reflection_table), 3))
    s0c[non_zero] = R.T @ s0
    # (R^T s1)^T = s1^T R for each row
    s1c[non_zero] = (s1[non_zero] / lengths[non_zero][:, np.newaxis]) @ R

    alignment_axis = (1.0, 0.0, 0.0)
    reflection_table["s0c"] = align_axis_along_z(
        alignment_axis, flumpy.vec_from_numpy(s0c)
    )
    reflection_table["s1c"] = align_axis_along_z(
        alignment_axis, flumpy.vec_from_numpy(s1c)
    )
    return reflection_table
