from __future__ import annotations

from dxtbx import flumpy
from dxtbx.model import tof_helpers

import dials_algorithms_spot_prediction_ext
from dials.array_family import flex
from dials_algorithms_spot_prediction_ext import (
    IndexGenerator,
    LaueRayPredictor,
    NaveStillsReflectionPredictor,
    PixelLabeller,
    PixelToMillerIndex,
    ReekeIndexGenerator,
    RotationAngles,
    ScanStaticRayPredictor,
    ScanVaryingRayPredictor,
    SphericalRelpStillsReflectionPredictor,
    StillsDeltaPsiReflectionPredictor,
    StillsRayPredictor,
    ray_intersection,
)

__all__ = [
    "IndexGenerator",
    "NaveStillsReflectionPredictor",
    "PixelLabeller",
    "PixelToMillerIndex",
    "ray_intersection",
    "ReekeIndexGenerator",
    "RotationAngles",
    "ScanStaticRayPredictor",
    "ScanStaticReflectionPredictor",
    "ScanVaryingRayPredictor",
    "ScanVaryingReflectionPredictor",
    "SphericalRelpStillsReflectionPredictor",
    "StillsDeltaPsiReflectionPredictor",
    "StillsRayPredictor",
    "StillsReflectionPredictor",
    "LaueRayPredictor",
    "LaueReflectionPredictor",
]


def ScanStaticReflectionPredictor(experiment, dmin=None, margin=1, padding=0, **kwargs):
    """
    A constructor for the reflection predictor.

    :param experiment: The experiment to predict for
    :param dmin: The maximum resolution to predict to
    :param margin: The margin for prediction
    :return: The spot predictor
    """

    # Get dmin if it is not set
    if dmin is None:
        dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

    # Only remove certain systematic absences
    space_group = experiment.crystal.get_space_group()
    space_group = space_group.build_derived_patterson_group()

    # Create the reflection predictor
    return dials_algorithms_spot_prediction_ext.ScanStaticReflectionPredictor(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        experiment.crystal.get_unit_cell(),
        space_group.type(),
        dmin,
        margin,
        padding,
    )


def ScanVaryingReflectionPredictor(
    experiment, dmin=None, margin=1, padding=0, **kwargs
):
    """
    A constructor for the reflection predictor.

    :param experiment: The experiment to predict for
    :param dmin: The maximum resolution to predict to
    :param margin: The margin for prediction
    :return: The spot predictor
    """

    # Get dmin if it is not set
    if dmin is None:
        dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

    # Only remove certain systematic absences
    space_group = experiment.crystal.get_space_group()
    space_group = space_group.build_derived_patterson_group()

    # Create the reflection predictor
    return dials_algorithms_spot_prediction_ext.ScanVaryingReflectionPredictor(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        space_group.type(),
        dmin,
        margin,
        padding,
    )


def StillsReflectionPredictor(experiment, dmin=None, spherical_relp=False, **kwargs):
    """
    A factory function for the reflection predictor.

    :param experiment: The experiment to predict for
    :param dmin: The maximum resolution to predict to
    :param spherical_relp: Whether to use the spherical relp prediction model
    :return: The spot predictor
    """

    # FIXME Selection of reflection predictor type is ugly. What is a better
    # way here? Should it be based entirely on the existence of certain types
    # of profile model within the experiment?

    # Get dmin if it is not set
    if dmin is None:
        dmin = experiment.detector.get_max_resolution(experiment.beam.get_s0())

    if spherical_relp:
        return SphericalRelpStillsReflectionPredictor(
            experiment.beam,
            experiment.detector,
            experiment.crystal.get_A(),
            experiment.crystal.get_unit_cell(),
            experiment.crystal.get_space_group().type(),
            dmin,
        )

    # Create the reflection predictor
    try:
        if (
            experiment.crystal.get_half_mosaicity_deg() is not None
            and experiment.crystal.get_domain_size_ang() is not None
        ):
            return NaveStillsReflectionPredictor(
                experiment.beam,
                experiment.detector,
                experiment.crystal.get_A(),
                experiment.crystal.get_unit_cell(),
                experiment.crystal.get_space_group().type(),
                dmin,
                experiment.crystal.get_half_mosaicity_deg(),
                experiment.crystal.get_domain_size_ang(),
            )
    except AttributeError:
        pass

    return StillsDeltaPsiReflectionPredictor(
        experiment.beam,
        experiment.detector,
        experiment.crystal.get_A(),
        experiment.crystal.get_unit_cell(),
        experiment.crystal.get_space_group().type(),
        dmin,
    )


def LaueReflectionPredictor(experiment, dmin: float):
    return dials_algorithms_spot_prediction_ext.LaueReflectionPredictor(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.crystal.get_A(),
        experiment.crystal.get_unit_cell(),
        experiment.crystal.get_space_group().type(),
        dmin,
    )


class TOFReflectionPredictor:
    def __init__(self, experiment, dmin):
        self.experiment = experiment
        self.dmin = dmin
        self.predictor = dials_algorithms_spot_prediction_ext.LaueReflectionPredictor(
            experiment.beam,
            experiment.detector,
            experiment.goniometer,
            experiment.crystal.get_A(),
            experiment.crystal.get_unit_cell(),
            experiment.crystal.get_space_group().type(),
            dmin,
        )

    def post_prediction(self, reflections):
        if "tof_cal" not in reflections:
            reflections["tof_cal"] = flex.double(reflections.nrows())
        if "L1" not in reflections:
            reflections["L1"] = flex.double(reflections.nrows())

        tof_cal = flex.double(reflections.nrows())
        L1 = flex.double(reflections.nrows())
        L0 = self.experiment.beam.get_sample_to_source_distance() * 10**-3  # (m)

        panel_numbers = flex.size_t(reflections["panel"])
        expt = self.experiment

        for i_panel in range(len(expt.detector)):
            sel = panel_numbers == i_panel
            expt_reflections = reflections.select(sel)
            x, y, _ = expt_reflections["xyzcal.mm"].parts()
            s1 = expt.detector[i_panel].get_lab_coord(flex.vec2_double(x, y))
            expt_L1 = s1.norms()
            expt_tof_cal = flex.double(expt_reflections.nrows())

            for idx in range(len(expt_reflections)):
                wavelength = expt_reflections[idx]["wavelength_cal"]
                tof = tof_helpers.tof_from_wavelength(
                    wavelength, L0 + expt_L1[idx] * 10**-3
                )
                expt_tof_cal[idx] = tof
            tof_cal.set_selected(sel, expt_tof_cal)
            L1.set_selected(sel, expt_L1)

        reflections["tof_cal"] = tof_cal
        reflections["L1"] = L1

        # Filter out predicted reflections outside of experiment range
        wavelength_range = expt.beam.get_wavelength_range()
        sel = reflections["wavelength_cal"] >= wavelength_range[0]
        reflections = reflections.select(sel)
        sel = reflections["wavelength_cal"] <= wavelength_range[1]
        reflections = reflections.select(sel)

        return reflections

    def indices_for_ub(self, indices):
        reflection_table = self.predictor(indices)
        reflection_table = self.post_prediction(reflection_table)

        interpolation_tof = self.experiment.scan.get_property("time_of_flight")
        interpolation_frames = list(range(len(interpolation_tof)))
        tof_to_frame = tof_helpers.tof_to_frame_interpolator(
            interpolation_tof, interpolation_frames
        )
        L0 = self.experiment.beam.get_sample_to_source_distance() * 10**-3  # (m)

        reflection_tof = (
            tof_helpers.tof_from_wavelength(
                reflection_table["wavelength_cal"],
                L0 + reflection_table["L1"] * 10**-3,
            )
            * 10**6
        )

        reflection_table = reflection_table.select(
            (reflection_tof > min(interpolation_tof))
            & (reflection_tof < max(interpolation_tof))
        )

        reflection_tof = reflection_tof.select(
            (reflection_tof > min(interpolation_tof))
            & (reflection_tof < max(interpolation_tof))
        )
        reflection_frames = flumpy.from_numpy(tof_to_frame(reflection_tof))
        x, y, _ = reflection_table["xyzcal.px"].parts()
        reflection_table["xyzcal.px"] = flex.vec3_double(x, y, reflection_frames)

        return reflection_table

    def for_ub(self, ub):
        reflection_table = self.predictor.for_ub(ub)
        reflection_table = self.post_prediction(reflection_table)

        interpolation_tof = self.experiment.scan.get_property("time_of_flight")
        interpolation_frames = list(range(len(interpolation_tof)))
        tof_to_frame = tof_helpers.tof_to_frame_interpolator(
            interpolation_tof, interpolation_frames
        )
        L0 = self.experiment.beam.get_sample_to_source_distance() * 10**-3  # (m)

        reflection_tof = (
            tof_helpers.tof_from_wavelength(
                reflection_table["wavelength_cal"],
                L0 + reflection_table["L1"] * 10**-3,
            )
            * 10**6
        )

        reflection_table = reflection_table.select(
            (reflection_tof > min(interpolation_tof))
            & (reflection_tof < max(interpolation_tof))
        )

        reflection_tof = reflection_tof.select(
            (reflection_tof > min(interpolation_tof))
            & (reflection_tof < max(interpolation_tof))
        )
        reflection_frames = flumpy.from_numpy(tof_to_frame(reflection_tof))
        x, y, _ = reflection_table["xyzcal.px"].parts()
        reflection_table["xyzcal.px"] = flex.vec3_double(x, y, reflection_frames)

        return reflection_table

    def for_reflection_table(self, reflections, UB):
        return self.predictor.for_reflection_table(reflections, UB)

    def all_reflections_for_asu(self, phi):
        return self.predictor.all_reflections_for_asu(float(phi))
