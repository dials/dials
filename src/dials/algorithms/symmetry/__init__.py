"""Methods for symmetry determination.

This module provides a base class for symmetry determination algorithms and
utility functions for symmetry analysis.
"""

from __future__ import annotations

import collections
import copy
import logging
import math
from io import StringIO

import libtbx
from cctbx import adptbx, sgtbx, uctbx
from cctbx.sgtbx.bravais_types import bravais_lattice
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from dxtbx.model import ExperimentList
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling, matthews

from dials.array_family import flex
from dials.util import resolution_analysis
from dials.util.exclude_images import (
    exclude_image_ranges_from_scans,
    get_selection_for_valid_image_ranges,
)
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.multi_dataset_handling import (
    select_datasets_on_identifiers,
    update_imageset_ids,
)
from dials.util.normalisation import quasi_normalisation

logger = logging.getLogger(__name__)


class symmetry_base:
    """Base class for symmetry analysis."""

    def __init__(
        self,
        intensities,
        normalisation="ml_aniso",
        lattice_symmetry_max_delta=2.0,
        d_min=libtbx.Auto,
        min_i_mean_over_sigma_mean=4,
        min_cc_half=0.6,
        relative_length_tolerance=None,
        absolute_angle_tolerance=None,
        best_monoclinic_beta=True,
        apply_sigma_correction=True,
    ):
        """Initialise a symmetry_base object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            symmetry analysis.
          normalisation (str): The normalisation method to use. Possible choices are
            'kernel', 'quasi', 'ml_iso' and 'ml_aniso'. Set to None to switch off
            normalisation altogether.
          lattice_symmetry_max_delta (float): The maximum value of delta for
            determining the lattice symmetry using the algorithm of Le Page (1982).
          d_min (float): Optional resolution cutoff to be applied to the input
            intensities. If set to :data:`libtbx.Auto` then d_min will be
            automatically determined according to the parameters
            ``min_i_mean_over_sigma_mean`` and ``min_cc_half``.
          min_i_mean_over_sigma_mean (float): minimum value of :math:`|I|/|sigma(I)|` for
            automatic determination of resolution cutoff.
          min_cc_half (float): minimum value of CC½ for automatic determination of
            resolution cutoff.
          relative_length_tolerance (float): Relative length tolerance in checking
            consistency of input unit cells against the median unit cell.
          absolute_angle_tolerance (float): Absolute angle tolerance in checking
            consistency of input unit cells against the median unit cell.
          best_monoclinic_beta (bool): If True, then for monoclinic centered cells, I2
            will be preferred over C2 if it gives a less oblique cell (i.e. smaller
            beta angle).
          apply_sigma_correction (bool): If True, correct SDs by "typical" SD factors.
        """
        self.input_intensities = intensities

        uc_params = [flex.double() for i in range(6)]
        for d in self.input_intensities:
            for i, p in enumerate(d.unit_cell().parameters()):
                uc_params[i].append(p)
        self.median_unit_cell = uctbx.unit_cell(
            parameters=[flex.median(p) for p in uc_params]
        )
        self._check_unit_cell_consistency(
            relative_length_tolerance, absolute_angle_tolerance
        )

        self.intensities = self.input_intensities[0]
        self.dataset_ids = flex.int(self.intensities.size(), 0)
        for i, d in enumerate(self.input_intensities[1:]):
            self.intensities = self.intensities.concatenate(
                d, assert_is_similar_symmetry=False
            )
            self.dataset_ids.extend(flex.int(d.size(), i + 1))
        self.intensities = self.intensities.customized_copy(
            unit_cell=self.median_unit_cell
        )
        self.intensities.set_observation_type_xray_intensity()
        sys_absent_flags = self.intensities.sys_absent_flags(integral_only=True).data()
        self.intensities = self.intensities.select(~sys_absent_flags)
        self.dataset_ids = self.dataset_ids.select(~sys_absent_flags)

        self.lattice_symmetry_max_delta = lattice_symmetry_max_delta
        self.subgroups = metric_subgroups(
            self.intensities.crystal_symmetry(),
            max_delta=self.lattice_symmetry_max_delta,
            bravais_types_only=False,
            best_monoclinic_beta=best_monoclinic_beta,
        )

        self.cb_op_inp_min = self.subgroups.cb_op_inp_minimum
        self.intensities = (
            self.intensities.change_basis(self.cb_op_inp_min)
            .customized_copy(space_group_info=sgtbx.space_group_info("P1"))
            .map_to_asu()
            .set_info(self.intensities.info())
        )

        self.lattice_group = (
            self.subgroups.result_groups[0]["subsym"].space_group().make_tidy()
        )
        self.patterson_group = (
            self.lattice_group.build_derived_patterson_group().make_tidy()
        )
        logger.info("Patterson group: %s", self.patterson_group.info())

        sel = self.patterson_group.epsilon(self.intensities.indices()) == 1
        self.intensities = self.intensities.select(sel)
        self.dataset_ids = self.dataset_ids.select(sel)

        # Correct SDs by "typical" SD factors
        if apply_sigma_correction:
            self._correct_sigmas(sd_fac=2.0, sd_b=0.0, sd_add=0.03)
        self._normalise(normalisation)
        self._resolution_filter(d_min, min_i_mean_over_sigma_mean, min_cc_half)

    def _check_unit_cell_consistency(
        self, relative_length_tolerance, absolute_angle_tolerance
    ):
        for d in self.input_intensities:
            if (
                relative_length_tolerance is not None
                and absolute_angle_tolerance is not None
            ):
                if not d.unit_cell().is_similar_to(
                    self.median_unit_cell,
                    relative_length_tolerance,
                    absolute_angle_tolerance,
                ):
                    raise ValueError(
                        f"Incompatible unit cell: {d.unit_cell()}\n"
                        + f"      median unit cell: {self.median_unit_cell}"
                    )

    def _normalise(self, method):
        if method is None:
            return
        elif method == "kernel":
            normalise = self.kernel_normalisation
        elif method == "quasi":
            normalise = quasi_normalisation
        elif method == "ml_iso":
            normalise = self.ml_iso_normalisation
        elif method == "ml_aniso":
            normalise = self.ml_aniso_normalisation

        normalised_intensities = None
        for i in range(int(flex.max(self.dataset_ids) + 1)):
            logger.info("\n" + "-" * 80 + "\n")
            logger.info("Normalising intensities for dataset %i\n", i + 1)
            intensities = self.intensities.select(self.dataset_ids == i)
            try:
                intensities = normalise(intensities)
            # Catch any of the several errors that can occur when there are too few
            # reflections for the selected normalisation routine.
            except (AttributeError, IndexError, RuntimeError):
                logger.warning(
                    "A problem occurred when trying to normalise the intensities by "
                    "the %s method.\n"
                    "There may be too few unique reflections. Resorting to using "
                    "un-normalised intensities instead.",
                    method,
                    exc_info=True,
                )
            if not normalised_intensities:
                normalised_intensities = intensities
            else:
                normalised_intensities = normalised_intensities.concatenate(intensities)
        self.intensities = normalised_intensities.set_info(
            self.intensities.info()
        ).set_observation_type_xray_intensity()

    def _correct_sigmas(self, sd_fac, sd_b, sd_add):
        # sd' = SDfac * Sqrt(sd^2 + SdB * I + (SDadd * I)^2)
        variance = flex.pow2(self.intensities.sigmas())
        si2 = flex.pow2(sd_add * self.intensities.data())
        ssc = variance + sd_b * self.intensities.data() + si2
        MINVARINFRAC = 0.1
        ssc.set_selected(ssc < MINVARINFRAC * variance, MINVARINFRAC * variance)
        sd = sd_fac * flex.sqrt(ssc)
        self.intensities = self.intensities.customized_copy(sigmas=sd).set_info(
            self.intensities.info()
        )

    @staticmethod
    def kernel_normalisation(intensities):
        """Kernel normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.
        """
        normalisation = absolute_scaling.kernel_normalisation(
            intensities, auto_kernel=True
        )
        return normalisation.normalised_miller.deep_copy().set_info(intensities.info())

    @staticmethod
    def ml_aniso_normalisation(intensities):
        """Anisotropic maximum-likelihood normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.
        """
        return symmetry_base._ml_normalisation(intensities, aniso=True)

    @staticmethod
    def ml_iso_normalisation(intensities):
        """Isotropic maximum-likelihood normalisation of the input intensities.

        Args:
          intensities (cctbx.miller.array): The intensities to be normalised.

        Returns:
          cctbx.miller.array: The normalised intensities.
        """
        return symmetry_base._ml_normalisation(intensities, aniso=False)

    @staticmethod
    def _ml_normalisation(intensities, aniso):
        # estimate number of residues per unit cell
        mr = matthews.matthews_rupp(intensities.crystal_symmetry())
        n_residues = mr.n_residues

        # estimate B-factor and scale factors for normalisation
        if aniso:
            normalisation = absolute_scaling.ml_aniso_absolute_scaling(
                intensities, n_residues=n_residues
            )
            if not normalisation.p_scale:
                raise RuntimeError("Unsuccessful normalisation")
            u_star = normalisation.u_star
        else:
            normalisation = absolute_scaling.ml_iso_absolute_scaling(
                intensities, n_residues=n_residues
            )
            if not (normalisation.b_wilson and normalisation.p_scale):
                raise RuntimeError("Unsuccessful normalisation")
            u_star = adptbx.b_as_u(
                adptbx.u_iso_as_u_star(intensities.unit_cell(), normalisation.b_wilson)
            )

        # record output in log file
        if aniso:
            b_cart = normalisation.b_cart
            logger.info("ML estimate of overall B_cart value:")
            logger.info(
                """\
  %5.2f, %5.2f, %5.2f
  %12.2f, %5.2f
  %19.2f""",
                b_cart[0],
                b_cart[3],
                b_cart[4],
                b_cart[1],
                b_cart[5],
                b_cart[2],
            )
        else:
            logger.info("ML estimate of overall B value:")
            logger.info("   %5.2f A**2", normalisation.b_wilson)
        logger.info("ML estimate of  -log of scale factor:")
        logger.info("  %5.2f", normalisation.p_scale)

        s = StringIO()
        mr.show(out=s)
        normalisation.show(out=s)
        logger.debug(s.getvalue())

        # apply scales
        return intensities.customized_copy(
            data=scaling.ml_normalise_aniso(
                intensities.indices(),
                intensities.data(),
                normalisation.p_scale,
                intensities.unit_cell(),
                u_star,
            ),
            sigmas=scaling.ml_normalise_aniso(
                intensities.indices(),
                intensities.sigmas(),
                normalisation.p_scale,
                intensities.unit_cell(),
                u_star,
            ),
        )

    def _resolution_filter(self, d_min, min_i_mean_over_sigma_mean, min_cc_half):
        logger.info("\n" + "-" * 80 + "\n")
        logger.info("Estimation of resolution for Laue group analysis\n")
        if d_min is libtbx.Auto and (
            min_i_mean_over_sigma_mean is not None or min_cc_half is not None
        ):
            d_min = resolution_filter_from_array(
                self.intensities, min_i_mean_over_sigma_mean, min_cc_half
            )
            if d_min is not None:
                logger.info("High resolution limit set to: %.2f", d_min)
            else:
                logger.info("High resolution limit set to: None")
        if d_min is not None:
            sel = self.intensities.resolution_filter_selection(d_min=d_min)
            self.intensities = self.intensities.select(sel).set_info(
                self.intensities.info()
            )
            self.dataset_ids = self.dataset_ids.select(sel)
            logger.info(
                "Selecting %i reflections with d > %.2f", self.intensities.size(), d_min
            )


def median_unit_cell(experiments):
    uc_params = [flex.double() for i in range(6)]
    for c in experiments.crystals():
        for i, p in enumerate(c.get_unit_cell().parameters()):
            uc_params[i].append(p)
    return uctbx.unit_cell(parameters=[flex.median(p) for p in uc_params])


def resolution_filter_from_array(intensities, min_i_mean_over_sigma_mean, min_cc_half):
    """Run the resolution filter using miller array data format."""
    rparams = resolution_analysis.phil_defaults.extract().resolution
    resolutionizer = resolution_analysis.Resolutionizer(intensities, rparams)
    return _resolution_filter(resolutionizer, min_i_mean_over_sigma_mean, min_cc_half)


def resolution_filter_from_reflections_experiments(
    reflections, experiments, min_i_mean_over_sigma_mean, min_cc_half
):
    """Run the resolution filter using native dials data formats."""
    rparams = resolution_analysis.phil_defaults.extract().resolution
    resolutionizer = (
        resolution_analysis.Resolutionizer.from_reflections_and_experiments(
            reflections, experiments, rparams
        )
    )
    return _resolution_filter(resolutionizer, min_i_mean_over_sigma_mean, min_cc_half)


def _resolution_filter(resolutionizer, min_i_mean_over_sigma_mean, min_cc_half):
    """Use a resolutionizer to perform filtering for symmetry analysis."""
    d_min_isigi = 0
    d_min_cc_half = 0
    if min_i_mean_over_sigma_mean is not None:
        try:
            d_min_isigi = resolutionizer.resolution(
                resolution_analysis.metrics.I_MEAN_OVER_SIGMA_MEAN,
                limit=min_i_mean_over_sigma_mean,
            ).d_min
        except RuntimeError as e:
            logger.info("I/σ(I) resolution filter failed with the following error:")
            logger.error(e)
        else:
            if d_min_isigi:
                logger.info(
                    "Resolution estimate from <I>/<σ(I)> > %.1f : %.2f",
                    min_i_mean_over_sigma_mean,
                    d_min_isigi,
                )
    if min_cc_half is not None:
        try:
            d_min_cc_half = resolutionizer.resolution(
                resolution_analysis.metrics.CC_HALF, limit=min_cc_half
            ).d_min
        except RuntimeError as e:
            logger.info("CC½ resolution filter failed with the following error:")
            logger.error(e)
        else:
            if d_min_cc_half:
                logger.info(
                    "Resolution estimate from CC½ > %.2f: %.2f",
                    min_cc_half,
                    d_min_cc_half,
                )
    valid = [d for d in (d_min_cc_half, d_min_isigi) if d]
    if valid:
        return min(valid)


def apply_change_of_basis_ops(experiments, reflections, change_of_basis_ops):
    """
    Apply the given change of basis ops to the input experiments and reflections

    Args:
        experiments (ExperimentList): a list of experiments.
        reflections (list): a list of reflection tables
        change_of_basis_ops (list): a list of cctbx.sgtbx.change_of_basis_op

    Returns: The experiments and reflections after application of the change of basis ops
    """

    for expt, refl, cb_op_inp_min in zip(experiments, reflections, change_of_basis_ops):
        refl["miller_index"] = cb_op_inp_min.apply(refl["miller_index"])
        expt.crystal = expt.crystal.change_basis(cb_op_inp_min)
        expt.crystal.set_space_group(sgtbx.space_group())
    return experiments, reflections


def eliminate_sys_absent(experiments, reflections):
    for i, expt in enumerate(experiments):
        if expt.crystal.get_space_group().n_ltr() > 1:
            effective_group = (
                expt.crystal.get_space_group().build_derived_reflection_intensity_group(
                    anomalous_flag=True
                )
            )
            sys_absent_flags = effective_group.is_sys_absent(
                reflections[i]["miller_index"]
            )
            if sys_absent_flags.count(True):
                reflections[i] = reflections[i].select(~sys_absent_flags)
                logger.info(
                    "Eliminating %i systematic absences for experiment %s",
                    sys_absent_flags.count(True),
                    expt.identifier,
                )
    return reflections


def get_subset_for_symmetry(experiments, reflection_tables, exclude_images=None):
    """Select an image range for symmetry analysis, or just select
    the first 360 degrees of data."""
    refls_for_sym = []
    if exclude_images:
        experiments = exclude_image_ranges_from_scans(
            reflection_tables, experiments, exclude_images
        )
        for refl, exp in zip(reflection_tables, experiments):
            sel = get_selection_for_valid_image_ranges(refl, exp)
            refls_for_sym.append(refl.select(sel))
    else:
        for expt, refl in zip(experiments, reflection_tables):
            sel = get_selection_for_valid_image_ranges(refl, expt)
            if expt.scan and not sel.count(False):
                # Use first 360 degrees if <360 deg i.e. first measured data,
                # but only if no reflections have been explicitly excluded
                # already
                scan_end = int(math.ceil(360 / abs(expt.scan.get_oscillation()[1])))
                if scan_end < len(expt.scan):
                    sel = refl["xyzobs.px.value"].parts()[2] <= scan_end
            refls_for_sym.append(refl.select(sel))
    return refls_for_sym


def unit_cells_are_similar_to(
    experiments, unit_cell, relative_length_tolerance, absolute_angle_tolerance
):
    return all(
        expt.crystal.get_unit_cell().is_similar_to(
            unit_cell,
            relative_length_tolerance=relative_length_tolerance,
            absolute_angle_tolerance=absolute_angle_tolerance,
        )
        for expt in experiments
    )


def change_of_basis_ops_to_minimum_cell(
    experiments, max_delta, relative_length_tolerance, absolute_angle_tolerance
):
    """
    Compute change of basis ops to map experiments to the minimum cell

    Map to the minimum cell via the best cell, which appears to guarantee that the
    resulting minimum cells are consistent.

    Args:
        experiments (ExperimentList): a list of experiments.
        reflections (list): a list of reflection tables

    Returns: The experiments and reflections mapped to the minimum cell
    """
    logger.info("Mapping all input cells to a common minimum cell")
    median_cell = median_unit_cell(experiments)
    unit_cells_are_similar = unit_cells_are_similar_to(
        experiments, median_cell, relative_length_tolerance, absolute_angle_tolerance
    )
    centring_symbols = [
        bravais_lattice(group=expt.crystal.get_space_group()).centring_symbol
        for expt in experiments
    ]
    if unit_cells_are_similar and len(set(centring_symbols)) == 1:
        groups = metric_subgroups(
            experiments[0]
            .crystal.get_crystal_symmetry()
            .customized_copy(unit_cell=median_cell),
            max_delta,
            enforce_max_delta_for_generated_two_folds=True,
        )
        group = groups.result_groups[0]
        cb_op_best_to_min = group["best_subsym"].change_of_basis_op_to_minimum_cell()
        cb_ops = [cb_op_best_to_min * group["cb_op_inp_best"]] * len(experiments)
    else:
        if len(set(centring_symbols)) > 1:
            logger.info(
                f"Multiple lattice centerings in input cells: {', '.join(list(set(centring_symbols)))}.\n"
                + "Attempting to map to a common minimum cell through the most common lattice group."
            )
        else:
            median_params = ", ".join(f"{i:.2f}" for i in median_cell.parameters())
            logger.info(
                f"Some input cells are not sufficiently similar to the median cell:\n  {median_params}\n"
                + f"  within a relative_length_tolerance of {relative_length_tolerance}\n"
                + f"  and an absolute_angle_tolerance of {absolute_angle_tolerance}.\n"
                + "Attempting to map to a common minimum cell through the most common lattice group."
            )
        groups = [
            metric_subgroups(
                expt.crystal.get_crystal_symmetry(),
                max_delta,
                best_monoclinic_beta=False,
                enforce_max_delta_for_generated_two_folds=True,
            )
            for expt in experiments
        ]
        counter = collections.Counter(
            g.result_groups[0]["best_subsym"].space_group() for g in groups
        )
        target_group = counter.most_common()[0][0]
        cb_ops = []
        best_cells = []
        for expt in experiments:
            groups = metric_subgroups(
                expt.crystal.get_crystal_symmetry(),
                max_delta,
                best_monoclinic_beta=False,
                enforce_max_delta_for_generated_two_folds=True,
            )
            group = None
            for g in groups.result_groups:
                if g["best_subsym"].space_group() == target_group:
                    group = g
            if group:
                cb_ops.append(group["cb_op_inp_best"])
                best_cells.append(group["best_subsym"].unit_cell())
            else:
                cb_ops.append(None)
                logger.info(
                    f"Couldn't match unit cell to target symmetry:\n"
                    f"{expt.crystal.get_crystal_symmetry()}\n"
                    f"Target symmetry: {target_group.info()}"
                )
        # now get median best cell
        from cctbx import uctbx

        uc_params = [flex.double() for i in range(6)]
        for unit_cell in best_cells:
            for i, p in enumerate(unit_cell.parameters()):
                uc_params[i].append(p)
        overall_best_unit_cell = uctbx.unit_cell(
            parameters=[flex.median(p) for p in uc_params]
        )
        n = 0
        for i, cb_op in enumerate(cb_ops):
            if cb_op is not None:
                best_cell = best_cells[n]
                n += 1
                if not best_cell.is_similar_to(
                    overall_best_unit_cell,
                    relative_length_tolerance,
                    absolute_angle_tolerance,
                ):
                    best_params = ", ".join(
                        f"{i:.2f}" for i in overall_best_unit_cell.parameters()
                    )
                    this_params = ", ".join(f"{i:.2f}" for i in best_cell.parameters())
                    logger.info(
                        f"Unable to map input cell for dataset {i} to a minimum cell through\n"
                        + f"a consistent best cell in space group {target_group.info()}.\n"
                        + f"Incompatible best cells:\n  {best_params},\n  {this_params},\n"
                        + f"  within a relative_length_tolerance of {relative_length_tolerance}\n"
                        + f"  and an absolute_angle_tolerance of {absolute_angle_tolerance}"
                    )
                    cb_ops[i] = None
        if not any(cb_ops):
            raise ValueError(
                "Exiting symmetry analysis: Unable to map any cells to a minimum cell through a consistent best cell"
            )

        ref_expts = ExperimentList(
            [expt for expt, cb_op in zip(experiments, cb_ops) if cb_op]
        ).change_basis(list(filter(None, cb_ops)))
        cb_op_ref_min = (
            ref_expts[0]
            .crystal.get_crystal_symmetry()
            .customized_copy(unit_cell=median_unit_cell(ref_expts))
            .change_of_basis_op_to_minimum_cell()
        )
        cb_ops = [cb_op_ref_min * cb_op if cb_op else None for cb_op in cb_ops]
    return cb_ops


def prepare_datasets_for_symmetry_analysis(
    experiments,
    reflection_tables,
    params,
):
    """Prepare datasets for symmetry analysis.

    This involves selecting suitable image ranges, mapping to a common minimum cell,
    eliminating systematic absences and filtering reflections by resolution.

    Args:
        experiments (ExperimentList): a list of experiments.
        reflection_tables (list): a list of reflection tables
        params (dials.command_line.cosym.phil_scope.extract()): The parameters for
            symmetry analysis.

    Returns:
        A tuple of (datasets,experiments, reflection_tables, cb_ops) after
        preparation for symmetry analysis.
    """
    # Map experiments and reflections to minimum cell
    cb_ops = change_of_basis_ops_to_minimum_cell(
        experiments,
        params.lattice_symmetry_max_delta,
        params.relative_length_tolerance,
        params.absolute_angle_tolerance,
    )
    exclude = [expt.identifier for expt, cb_op in zip(experiments, cb_ops) if not cb_op]
    if len(exclude):
        if not params.exclude_inconsistent_unit_cells:
            indices = [str(i + 1) for i, v in enumerate(cb_ops) if v is None]
            raise ValueError(
                "Exiting symmetry analysis: Unable to match all cells to target symmetry.\n"
                + "This may be avoidable by increasing the absolute_angle_tolerance or relative_length_tolerance,\n"
                + "if the cells are similar enough.\n"
                + f"Alternatively, remove dataset number{'s' if len(indices) > 1 else ''} {', '.join(indices)} from the input"
            )
        exclude_indices = [i for i, cb_op in enumerate(cb_ops) if not cb_op]
        logger.info(
            f"Excluding {len(exclude)} datasets from further analysis "
            f"(couldn't determine consistent cb_op to minimum cell):\n"
            f"dataset indices: {exclude_indices}",
        )
        logger.info(
            "This may be avoidable by increasing the absolute_angle_tolerance or relative_length_tolerance,\n"
            + "if the cells are similar enough.\n"
        )

        if params.output.excluded:
            excluded_experiments = copy.deepcopy(experiments)
            excluded_reflections = copy.deepcopy(reflection_tables)
            excluded_experiments, excluded_reflections = select_datasets_on_identifiers(
                excluded_experiments, excluded_reflections, use_datasets=exclude
            )
            logger.info(
                "Saving excluded experiments to %s",
                params.output.excluded_prefix + ".expt",
            )
            excluded_experiments.as_file(params.output.excluded_prefix + ".expt")
            logger.info(
                "Saving excluded reflections to %s",
                params.output.excluded_prefix + ".refl",
            )
            joined = flex.reflection_table()
            excluded_reflections = update_imageset_ids(
                excluded_experiments, excluded_reflections
            )
            for refl in excluded_reflections:
                joined.extend(refl)
            joined.as_file(params.output.excluded_prefix + ".refl")

        experiments, reflection_tables = select_datasets_on_identifiers(
            experiments, reflection_tables, exclude_datasets=exclude
        )
        cb_ops = list(filter(None, cb_ops))

    # Eliminate reflections that are systematically absent due to centring
    # of the lattice, otherwise they would lead to non-integer miller indices
    # when reindexing to a primitive setting
    reflection_tables = eliminate_sys_absent(experiments, reflection_tables)

    experiments, reflection_tables = apply_change_of_basis_ops(
        experiments, reflection_tables, cb_ops
    )

    refls_for_sym = get_subset_for_symmetry(
        experiments, reflection_tables, params.exclude_images
    )

    # Transform models into miller arrays
    datasets = filtered_arrays_from_experiments_reflections(
        experiments,
        refls_for_sym,
        outlier_rejection_after_filter=True,
        partiality_threshold=params.partiality_threshold,
    )

    datasets = [ma.as_anomalous_array().merge_equivalents().array() for ma in datasets]

    return datasets, experiments, refls_for_sym, cb_ops
