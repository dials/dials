"""
Collection of factories for creating the scalers.
"""

from __future__ import annotations

import logging

from libtbx import Auto

from dials.algorithms.scaling.scaler import (
    MultiScaler,
    NullScaler,
    SingleScaler,
    TargetScaler,
)
from dials.algorithms.scaling.scaling_library import choose_initial_scaling_intensities
from dials.algorithms.scaling.scaling_utilities import (
    BadDatasetForScalingException,
    Reasons,
    align_axis_along_z,
    calc_crystal_frame_vectors,
    quasi_normalisation,
)
from dials.array_family import flex
from dials.util.filter_reflections import (
    filter_reflection_table_selection,
    sum_partial_reflections,
)

logger = logging.getLogger("dials")


def create_scaler(params, experiments, reflections):
    """Read an experimentlist and list of reflection tables and return
    an appropriate scaler. Requires experiment identifiers are correctly set in
    the experiments and reflections."""
    if not reflections:
        raise ValueError("No reflection tables provided as input")
    if len(reflections) == 1:
        scaler = SingleScalerFactory.create(params, experiments[0], reflections[0])
    else:
        is_scaled_list = [expt.scaling_model.is_scaled for expt in experiments]
        # if target mtz/model -> want to do targeted scaling only
        if params.scaling_options.target_mtz or params.scaling_options.target_model:
            # last experiment/refl is target, rest are to scale against this
            scaler = TargetScalerFactory.create_for_target_against_reference(
                params, experiments, reflections
            )
        elif (  # not mtz/model, but some are scaled
            params.scaling_options.target_cycle and len(set(is_scaled_list)) == 2
        ):  # if only some scaled and want to do targeted scaling
            scaler = TargetScalerFactory.create(params, experiments, reflections)
        else:  # else just make one multiscaler for all refls
            scaler = MultiScalerFactory.create(params, experiments, reflections)
    return scaler


class ScalerFactory:
    """Base class for Scaler Factories"""

    @staticmethod
    def filter_bad_reflections(
        reflections, partiality_cutoff=0.4, min_isigi=-5.0, intensity_choice="combine"
    ):
        """Initial filter to select integrated reflections."""
        logger.info(
            "Applying filter of min_isigi > %s, partiality > %s",
            min_isigi,
            partiality_cutoff,
        )
        logger.disabled = True
        if intensity_choice == "combine":
            if "intensity.sum.value" not in reflections:
                if "intensity.prf.value" not in reflections:
                    intensity_choice = None
                else:
                    intensity_choice = ["profile"]
            elif "intensity.prf.value" not in reflections:
                intensity_choice = ["sum"]
            else:
                intensity_choice = ["sum | profile"]
        else:
            intensity_choice = [intensity_choice]
        if intensity_choice:
            good = filter_reflection_table_selection(
                reflections,
                intensity_choice=intensity_choice,
                combine_partials=False,
                partiality_threshold=partiality_cutoff,
                min_isigi=min_isigi,
            )
            mask = ~good
            reflections.set_flags(mask, reflections.flags.excluded_for_scaling)
        logger.disabled = False
        return reflections

    @staticmethod
    def ensure_experiment_identifier(experiment, reflection_table):
        """Check for consistent experiment identifier, and if not set then set it
        using scaled_id."""
        id_vals = list(reflection_table.experiment_identifiers().values())
        assert experiment.identifier in id_vals, (experiment.identifier, list(id_vals))
        assert len(id_vals) == 1, list(id_vals)
        logger.info(
            "The experiment id for this dataset is %s.",
            reflection_table.experiment_identifiers().keys()[0],
        )


class SingleScalerFactory(ScalerFactory):
    """Factory for creating a scaler for a single dataset."""

    @classmethod
    def create(cls, params, experiment, reflection_table, for_multi=False):
        """Perform reflection_table preprocessing and create a SingleScaler."""

        cls.ensure_experiment_identifier(experiment, reflection_table)

        logger.info(
            "The scaling model type being applied is %s. \n",
            experiment.scaling_model.id_,
        )
        try:
            reflection_table = cls.filter_bad_reflections(
                reflection_table,
                partiality_cutoff=params.cut_data.partiality_cutoff,
                min_isigi=params.cut_data.min_isigi,
                intensity_choice=params.reflection_selection.intensity_choice,
            )
        except ValueError:
            raise BadDatasetForScalingException

        # combine partial measurements of same reflection, to handle those reflections
        # that were split by dials.integrate  - changes size of reflection table.
        reflection_table = sum_partial_reflections(reflection_table)

        if "inverse_scale_factor" not in reflection_table:
            reflection_table["inverse_scale_factor"] = flex.double(
                reflection_table.size(), 1.0
            )
        elif (
            reflection_table["inverse_scale_factor"].count(0.0)
            == reflection_table.size()
        ):
            reflection_table["inverse_scale_factor"] = flex.double(
                reflection_table.size(), 1.0
            )
        reflection_table = choose_initial_scaling_intensities(
            reflection_table, params.reflection_selection.intensity_choice
        )

        excluded_for_scaling = reflection_table.get_flags(
            reflection_table.flags.excluded_for_scaling
        )
        user_excluded = reflection_table.get_flags(
            reflection_table.flags.user_excluded_in_scaling
        )
        reasons = Reasons()
        reasons.add_reason("user excluded", user_excluded.count(True))
        reasons.add_reason("excluded for scaling", excluded_for_scaling.count(True))
        n_excluded = (excluded_for_scaling | user_excluded).count(True)
        if n_excluded == reflection_table.size():
            logger.info("All reflections were determined to be unsuitable for scaling.")
            logger.info(reasons)
            raise BadDatasetForScalingException(
                """Unable to use this dataset for scaling"""
            )
        else:
            logger.info(
                "Excluding %s/%s reflections\n%s",
                n_excluded,
                reflection_table.size(),
                reasons,
            )

        if params.reflection_selection.method == "intensity_ranges":
            reflection_table = quasi_normalisation(reflection_table, experiment)
        if (
            params.reflection_selection.method in (None, Auto, "auto", "quasi_random")
        ) or (
            experiment.scaling_model.id_ == "physical"
            and "absorption" in experiment.scaling_model.components
        ):
            if experiment.scan:
                reflection_table = calc_crystal_frame_vectors(
                    reflection_table, experiment
                )
                alignment_axis = (1.0, 0.0, 0.0)
                reflection_table["s0c"] = align_axis_along_z(
                    alignment_axis, reflection_table["s0c"]
                )
                reflection_table["s1c"] = align_axis_along_z(
                    alignment_axis, reflection_table["s1c"]
                )
        try:
            scaler = SingleScaler(params, experiment, reflection_table, for_multi)
        except BadDatasetForScalingException as e:
            raise ValueError(e)
        else:
            return scaler


class NullScalerFactory(ScalerFactory):
    "Factory for creating null scaler"

    @classmethod
    def create(cls, params, experiment, reflection_table):
        """Return Null Scaler."""

        logger.info("Preprocessing target dataset for scaling. \n")
        reflection_table = cls.filter_bad_reflections(reflection_table)
        if "variance" in reflection_table:
            variance_mask = reflection_table["variance"] <= 0.0
            reflection_table.set_flags(
                variance_mask, reflection_table.flags.excluded_for_scaling
            )
        else:
            reflection_table.set_flags(
                flex.bool(reflection_table.size(), False),
                reflection_table.flags.excluded_for_scaling,
            )
        logger.info(
            "%s reflections not suitable for scaling\n",
            reflection_table.get_flags(
                reflection_table.flags.excluded_for_scaling
            ).count(True),
        )
        cls.ensure_experiment_identifier(experiment, reflection_table)
        return NullScaler(params, experiment, reflection_table)


class MultiScalerFactory:
    "Factory for creating a scaler for multiple datasets"

    @staticmethod
    def create(params, experiments, reflections):
        """create a list of single scalers to pass to a MultiScaler."""
        single_scalers = []
        idx_to_remove = []
        for i, (expt, refl) in enumerate(zip(experiments, reflections)):
            # Remove bad datasets that literally have no integrated reflections
            try:
                scaler = SingleScalerFactory.create(params, expt, refl, for_multi=True)
            except BadDatasetForScalingException as e:
                logger.info(e)
                idx_to_remove.append(i)
            else:
                single_scalers.append(scaler)
        ms = MultiScaler(single_scalers)
        if idx_to_remove:
            for i in idx_to_remove:
                ms.removed_datasets.append(experiments[i].identifier)
            logger.info(
                "Removed experiments %s", " ".join(str(i) for i in idx_to_remove)
            )
        return ms

    @staticmethod
    def create_from_targetscaler(targetscaler):
        """method to pass scalers from TargetScaler to a MultiScaler"""
        single_scalers = []
        for scaler in targetscaler.unscaled_scalers:
            single_scalers.append(scaler)
        single_scalers.extend(targetscaler.single_scalers)
        multiscaler = MultiScaler(single_scalers)
        multiscaler.observers = targetscaler.observers
        return multiscaler


class TargetScalerFactory:
    "Factory for creating a targeted scaler for multiple datasets"

    @staticmethod
    def create_for_target_against_reference(params, experiments, reflections):
        """Create TargetScaler for case where have a target_mtz or target_model."""
        scaled_scalers = []
        unscaled_scalers = []
        idx_to_remove = []

        for i, (expt, refl) in enumerate(zip(experiments[:-1], reflections[:-1])):
            # Remove bad datasets that literally have no integrated reflections
            try:
                scaler = SingleScalerFactory.create(params, expt, refl, for_multi=True)
            except BadDatasetForScalingException as e:
                logger.info(e)
                idx_to_remove.append(i)
            else:
                unscaled_scalers.append(scaler)
        scaled_scalers = [
            NullScalerFactory.create(params, experiments[-1], reflections[-1])
        ]
        ts = TargetScaler(scaled_scalers, unscaled_scalers)
        if idx_to_remove:
            for i in idx_to_remove:
                ts.removed_datasets.append(experiments[i].identifier)
            logger.info(
                "Removed experiments %s", " ".join(str(i) for i in idx_to_remove)
            )
        return ts

    @staticmethod
    def create(params, experiments, reflections):
        """sort scaled and unscaled datasets to pass to TargetScaler"""
        scaled_experiments = []
        scaled_scalers = []
        unscaled_scalers = []
        idx_to_remove = []

        for i, (expt, refl) in enumerate(zip(experiments, reflections)):
            # Remove bad datasets that literally have no integrated reflections
            try:
                scaler = SingleScalerFactory.create(params, expt, refl, for_multi=True)
            except BadDatasetForScalingException as e:
                logger.info(e)
                idx_to_remove.append(i)
            else:
                if expt.scaling_model.is_scaled:
                    scaled_scalers.append(scaler)
                    scaled_experiments.append(expt)
                else:
                    unscaled_scalers.append(scaler)
        if idx_to_remove:
            for j in idx_to_remove[::-1]:
                del experiments[j]
                del reflections[j]
            logger.info(
                "Removed experiments %s", " ".join(str(i) for i in idx_to_remove)
            )

        n_exp, n_refl = (len(experiments), len(reflections))
        n_ss, n_us = (len(scaled_scalers), len(unscaled_scalers))
        assert n_exp == n_ss + n_us, (n_exp, str(n_ss) + " + " + str(n_us))
        assert n_exp == n_refl, (n_exp, n_refl)
        return TargetScaler(scaled_scalers, unscaled_scalers)
