"""
This module defines classes which implement the stages of the scaling algorithm.

These 'scalers' act to initialise and connect various parts of the scaling
algorithm and datastructures such as the Ih_table etc, and
present a united interface to the main scaling algorithm for single, multi
and targeted scaling.

The SingleScaler is defined, for scaling of a single dataset, a MultiScaler is
defined for scaling multiple datasets simultaneously and a TargetScaler is
defined for targeted scaling.
"""

from __future__ import annotations

import copy
import logging
import time
from io import StringIO
from math import ceil

import numpy as np

from dxtbx import flumpy
from libtbx import Auto
from scitbx import sparse

from dials.algorithms.scaling.basis_functions import RefinerCalculator
from dials.algorithms.scaling.combine_intensities import (
    MultiDatasetIntensityCombiner,
    SingleDatasetIntensityCombiner,
)
from dials.algorithms.scaling.error_model.engine import run_error_model_refinement
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.outlier_rejection import (
    determine_Esq_outlier_index_arrays,
    determine_outlier_index_arrays,
)
from dials.algorithms.scaling.parameter_handler import ScalingParameterManagerGenerator
from dials.algorithms.scaling.reflection_selection import (
    _select_groups_on_Isigma_cutoff,
    calculate_scaling_subset_ranges,
    calculate_scaling_subset_ranges_with_E2,
    select_connected_reflections_across_datasets,
)
from dials.algorithms.scaling.scaling_library import (
    merging_stats_from_scaled_array,
    scaled_data_as_miller_array,
)
from dials.algorithms.scaling.scaling_refiner import scaling_refinery
from dials.algorithms.scaling.scaling_utilities import (
    BadDatasetForScalingException,
    DialsMergingStatisticsError,
    log_memory_usage,
)
from dials.algorithms.scaling.target_function import ScalingTarget, ScalingTargetFixedIH
from dials.array_family import flex
from dials.util import tabulate
from dials.util.observer import Subject
from dials_scaling_ext import calc_sigmasq as cpp_calc_sigmasq
from dials_scaling_ext import row_multiply

logger = logging.getLogger("dials")

flex.set_random_seed(42)


class ScalerBase(Subject):
    """
    Abstract base class for all scalers (single and multiple).
    """

    def __init__(self, params):
        """Define the properties of a scaler."""
        super().__init__(
            events=[
                "performed_scaling",
                "performed_error_analysis",
                "performed_outlier_rejection",
            ]
        )
        self._params = params
        self._Ih_table = None
        self._global_Ih_table = None
        self._free_Ih_table = None
        self._work_free_stats = []
        self._removed_datasets = []
        self._active_scalers = []

    @property
    def active_scalers(self):
        """A list of scalers that are currently being used in the algorithm."""
        return self._active_scalers

    @property
    def removed_datasets(self):
        """The list of removed datasets."""
        return self._removed_datasets

    @property
    def work_free_stats(self):
        """Holder for work/free set statistics."""
        return self._work_free_stats

    @property
    def Ih_table(self):
        """The Ih_table datastructure for use in minimisation."""
        return self._Ih_table

    @property
    def global_Ih_table(self):
        """
        An Ih_table datastructure containing all suitable reflections.

        This includes reflections across all datasets being minimised, and there
        should only be one instance, maintained by the highest level scaler, e.g.
        a multiscaler in a multi-dataset case.
        """
        return self._global_Ih_table

    @property
    def params(self):
        """The params phil scope."""
        return self._params

    ### Interface for scaling refiner

    def update_for_minimisation(self, apm, block_id):
        """Update the scale factors and Ih for the next minimisation iteration."""
        raise NotImplementedError()

    def get_blocks_for_minimisation(self):
        """Return the blocks to iterate over during refinement."""
        if self.Ih_table.free_Ih_table:
            return self.Ih_table.blocked_data_list[:-1]
        return self.Ih_table.blocked_data_list

    ## Interface for algorithms using the scaler

    def round_of_outlier_rejection(self):
        """Perform a round of outlier rejection on the reflections."""
        raise NotImplementedError()

    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
        """Expand scales from a subset to all reflections."""
        raise NotImplementedError()

    def combine_intensities(self):
        """Combine prf and sum intensities to give optimal intensities."""
        pass

    def make_ready_for_scaling(self, outlier=True):
        """Make the scaler in a prepared state for scaling."""
        raise NotImplementedError()

    def determine_shared_model_components(self):
        return None

    @Subject.notify_event(event="performed_scaling")
    def perform_scaling(self, engine=None, max_iterations=None, tolerance=None):
        """Minimise the scaling model."""
        self._perform_scaling(
            target_type=ScalingTarget,
            engine=engine,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )

    def _perform_scaling(
        self, target_type, engine=None, max_iterations=None, tolerance=None
    ):
        target = target_type()
        # find if any components to share
        shared = self.determine_shared_model_components()
        pmg = ScalingParameterManagerGenerator(
            self.active_scalers,
            target,
            self.params.scaling_refinery.refinement_order,
            shared=shared,
        )
        for apm in pmg.parameter_managers():
            if not engine:
                engine = self.params.scaling_refinery.engine
            if not max_iterations:
                max_iterations = self.params.scaling_refinery.max_iterations
            st = time.time()
            refinery = scaling_refinery(
                engine=engine,
                scaler=self,
                target=target,
                prediction_parameterisation=apm,
                max_iterations=max_iterations,
            )
            if tolerance:
                refinery.set_tolerance(tolerance)
            try:
                refinery.run()
            except RuntimeError as e:
                logger.error(e, exc_info=True)
            logger.info("Time taken for refinement %.2f", (time.time() - st))
            refinery.print_step_table()
            self._update_after_minimisation(apm)
            logger.info("\n" + "=" * 80 + "\n")

    def clear_Ih_table(self):
        """Delete the data from the current Ih_table."""
        self._Ih_table = []

    def fix_initial_parameter(self):
        return False

    def prepare_reflection_tables_for_output(self):
        """Finish adjust reflection table data at the end of the algorithm."""
        # First adjust variances
        for scaler in self.active_scalers:
            if scaler._experiment.scaling_model.error_model:
                scaler.reflection_table["variance"] = flumpy.from_numpy(
                    scaler._experiment.scaling_model.error_model.update_variances(
                        scaler.reflection_table["variance"],
                        scaler.reflection_table["intensity"],
                    )
                )
            # now increase the errors slightly to take into account the uncertainty in the
            # inverse scale factors
            if (
                scaler.var_cov_matrix.non_zeroes > 0
            ):  # parameters errors have been determined
                fractional_error = (
                    flex.sqrt(scaler.reflection_table["inverse_scale_factor_variance"])
                    / scaler.reflection_table["inverse_scale_factor"]
                )
                variance_scaling = (
                    flex.double(scaler.reflection_table.size(), 1.0) + fractional_error
                )
                scaler.reflection_table["variance"] *= variance_scaling
        self._set_outliers()
        self._clean_reflection_tables()
        msg = """
The reflection table variances have been adjusted to account for the
uncertainty in the scaling model""" + (
            "s for all datasets" if len(self.active_scalers) > 1 else ""
        )
        logger.info(msg)
        if self._free_Ih_table:
            # calc merging stats and log.
            free_miller_array = scaled_data_as_miller_array(
                [self.get_free_set_reflections()],
                [self.active_scalers[0].experiment],
                anomalous_flag=False,
            )
            work_miller_array = scaled_data_as_miller_array(
                [self.get_work_set_reflections()],
                [self.active_scalers[0].experiment],
                anomalous_flag=False,
            )
            free, _ = merging_stats_from_scaled_array(
                free_miller_array,
                self.params.output.merging.nbins,
                self.params.output.use_internal_variance,
                anomalous=False,
            )
            s = StringIO()
            s.write("\nFree set statistics\n")
            free.show(out=s)
            work, _ = merging_stats_from_scaled_array(
                work_miller_array,
                self.params.output.merging.nbins,
                self.params.output.use_internal_variance,
                anomalous=False,
            )
            s.write("\nWork set statistics\n")
            work.show(out=s)
            logger.debug(s.getvalue())
            self._work_free_stats = [
                work.overall.r_meas,
                free.overall.r_meas,
                free.overall.r_meas - work.overall.r_meas,
                work.overall.cc_one_half,
                free.overall.cc_one_half,
                work.overall.cc_one_half - free.overall.cc_one_half,
            ]

    # Internal general scaler methods

    def _set_outliers(self):
        """Set the scaling outliers in the individual reflection tables."""
        for scaler in self.active_scalers:
            suitable_isel = scaler.suitable_refl_for_scaling_sel.iselection()
            outlier_isel = suitable_isel.select(scaler.outliers)
            n = scaler.reflection_table.size()
            outliers_mask = flex.bool(n, False)
            outliers_mask.set_selected(outlier_isel, True)
            scaler.reflection_table.unset_flags(
                flex.bool(n, True), scaler.reflection_table.flags.outlier_in_scaling
            )
            scaler.reflection_table.set_flags(
                outliers_mask, scaler.reflection_table.flags.outlier_in_scaling
            )

    def _clean_reflection_tables(self):
        """Remove unnecessary columns added to reflection tables."""
        for scaler in self.active_scalers:
            scaler.clean_reflection_table()

    def _update_after_minimisation(self, parameter_manager):
        if parameter_manager.apm_list[0].var_cov_matrix:
            for i, scaler in enumerate(self.active_scalers):
                scaler.update_var_cov(parameter_manager.apm_list[i])


class SingleScaler(ScalerBase):
    """Definition of a scaler for a single dataset."""

    id_ = "single"

    def __init__(self, params, experiment, reflection_table, for_multi=False):
        """
        Initialise a single-dataset scaler.

        The reflection table needs the columns 'inverse_scale_factor', 'Esq',
        'intensity', 'variance', 'id', which are guaranteed if the scaler is
        created using the SingleScalerFactory.
        """
        assert all(
            i in reflection_table
            for i in ["inverse_scale_factor", "intensity", "variance", "id"]
        )
        super().__init__(params)
        self._experiment = experiment
        n_model_params = sum(val.n_params for val in self.components.values())
        self._var_cov_matrix = sparse.matrix(n_model_params, n_model_params)
        self._initial_keys = list(reflection_table.keys())
        self._reflection_table = reflection_table
        self._Ih_table = None  # stores data for reflections used for minimisation
        self.suitable_refl_for_scaling_sel = self._get_suitable_for_scaling_sel(
            self._reflection_table
        )
        self.n_suitable_refl = self.suitable_refl_for_scaling_sel.count(True)
        if self._experiment.scaling_model.is_scaled:
            outliers = self._reflection_table.get_flags(
                self._reflection_table.flags.outlier_in_scaling
            )
            self.outliers = outliers.select(self.suitable_refl_for_scaling_sel)
        else:
            self.outliers = flex.bool(self.n_suitable_refl, False)
        self.scaling_subset_sel = (
            None  # A selection of len n_suitable_refl of scaling subset selection
        )
        self.scaling_selection = None  # As above, but with outliers deselected also
        self.free_set_selection = flex.bool(self.n_suitable_refl, False)
        self._free_Ih_table = None  # An array of len n_suitable_refl
        self._configure_model_and_datastructures(for_multi=for_multi)
        if self.params.weighting.error_model.error_model:
            # reload current error model parameters, or create new null
            self.experiment.scaling_model.load_error_model(
                self.params.weighting.error_model
            )
            self._update_error_model(self.experiment.scaling_model.error_model)
        if "Imid" in self.experiment.scaling_model.configdict:
            self._combine_intensities(self.experiment.scaling_model.configdict["Imid"])
        if not self._experiment.scaling_model.is_scaled:
            self.round_of_outlier_rejection()
        if not for_multi:
            self._select_reflections_for_scaling()
            self._create_Ih_table()
            self._update_model_data()
        else:
            self._global_Ih_table = None
            self.scaling_selection = ~self.outliers
        logger.info(
            "Completed preprocessing and initialisation for this dataset.\n"
            "\n" + "=" * 80 + "\n"
        )
        log_memory_usage()

    @Subject.notify_event(event="performed_error_analysis")
    def perform_error_optimisation(self, update_Ih=True):
        """Perform an optimisation of the sigma values."""
        # error model should be determined using anomalous groups
        Ih_table, _ = self._create_global_Ih_table(anomalous=True, remove_outliers=True)
        try:
            model = run_error_model_refinement(
                self._experiment.scaling_model.error_model, Ih_table
            )
        except (ValueError, RuntimeError) as e:
            logger.info(e)
            logger.debug(e, exc_info=True)
        else:
            self._update_error_model(model, update_Ih=update_Ih)

    def _update_error_model(self, error_model, update_Ih=True):
        """Update the error model in Ih table."""
        if update_Ih and error_model:
            self.global_Ih_table.update_weights(error_model)
            if self._free_Ih_table:
                self._free_Ih_table.update_weights(error_model)
        self._experiment.scaling_model.set_error_model(error_model)

    def get_valid_reflections(self):
        """All reflections not bad for scaling or user excluded."""
        return self.reflection_table.select(self.suitable_refl_for_scaling_sel)

    def get_work_set_reflections(self):
        valid = self.get_valid_reflections()
        return valid.select(~self.free_set_selection)

    def get_free_set_reflections(self):
        """Get all reflections in the free set if it exists."""
        valid = self.get_valid_reflections()
        return valid.select(self.free_set_selection)

    def get_reflections_for_model_minimisation(self):
        valid = self.get_valid_reflections()
        return valid.select(self.scaling_selection)

    @property
    def active_scalers(self):
        return [self]

    @property
    def experiment(self):
        """The experiment object associated with the dataset."""
        return self._experiment

    @property
    def reflection_table(self):
        """The reflection table of the dataset."""
        return self._reflection_table

    @reflection_table.setter
    def reflection_table(self, new_table):
        """Set the reflection table of the dataset."""
        self._reflection_table = new_table

    def fix_initial_parameter(self):
        fixed = self.experiment.scaling_model.fix_initial_parameter(self.params)
        return fixed

    @staticmethod
    def _get_suitable_for_scaling_sel(reflections):
        """Extract suitable reflections for scaling from the reflection table."""
        user_excl = reflections.get_flags(reflections.flags.user_excluded_in_scaling)
        excl_for_scale = reflections.get_flags(reflections.flags.excluded_for_scaling)
        suitable_refl_for_scaling_sel = ~(user_excl | excl_for_scale)
        return suitable_refl_for_scaling_sel

    @property
    def components(self):
        """Shortcut to the scaling model components."""
        return self.experiment.scaling_model.components

    @property
    def fixed_components(self):
        """Shortcut to the scaling model fixed components."""
        return self.experiment.scaling_model.fixed_components

    @property
    def consecutive_refinement_order(self):
        """Link to consecutive refinement order for parameter manager."""
        return self.experiment.scaling_model.consecutive_refinement_order

    @property
    def var_cov_matrix(self):
        """The variance covariance matrix for the parameters."""
        return self._var_cov_matrix

    def update_var_cov(self, apm):
        """
        Update the full parameter variance covariance matrix after a refinement.

        If all parameters have been refined, then the full var_cov matrix can be set.
        Else one must select subblocks for pairs of parameters and assign these into
        the full var_cov matrix, taking care to out these in the correct position.
        This is applicable if only some parameters have been refined in this cycle.
        """
        var_cov_list = apm.var_cov_matrix  # values are passed as a list from refinery
        if int(var_cov_list.size() ** 0.5) == self.var_cov_matrix.n_rows:
            self._var_cov_matrix.assign_block(
                var_cov_list.matrix_copy_block(
                    0, 0, apm.n_active_params, apm.n_active_params
                ),
                0,
                0,
            )
        else:  # need to set part of the var_cov matrix e.g. if only refined some params
            # first work out the order in self._var_cov_matrix
            cumul_pos_dict = {}
            n_cumul_params = 0
            for name, component in self.components.items():
                cumul_pos_dict[name] = n_cumul_params
                n_cumul_params += component.n_params
            # now get a var_cov_matrix subblock for pairs of parameters
            for name in apm.components_list:
                for name2 in apm.components_list:
                    n_rows = apm.components[name]["n_params"]
                    n_cols = apm.components[name2]["n_params"]
                    start_row = apm.components[name]["start_idx"]
                    start_col = apm.components[name2]["start_idx"]
                    sub = var_cov_list.matrix_copy_block(
                        start_row, start_col, n_rows, n_cols
                    )
                    # now set this block into correct location in overall var_cov
                    self._var_cov_matrix.assign_block(
                        sub, cumul_pos_dict[name], cumul_pos_dict[name2]
                    )

    def update_for_minimisation(self, apm, block_id):
        """Update the scale factors and Ih for the next minimisation iteration."""
        apm_i = apm.apm_list[0]
        scales_i, derivs_i = RefinerCalculator.calculate_scales_and_derivatives(
            apm_i, block_id
        )
        self.Ih_table.set_derivatives(derivs_i, block_id)
        self.Ih_table.set_inverse_scale_factors(flumpy.to_numpy(scales_i), block_id)
        self.Ih_table.calc_Ih(block_id)

    def combine_intensities(self):
        """Combine prf and sum intensities to give optimal intensities."""
        self._combine_intensities()

    def _combine_intensities(self, use_Imid=None):
        try:
            if use_Imid is not None:
                logger.info(
                    "Using previously determined optimal intensity choice: %s\n",
                    {
                        use_Imid: str(round(use_Imid, 4)),
                        0: "profile intensities",
                        1: "summation intensities",
                    }[use_Imid],
                )
            else:
                logger.info("Performing profile/summation intensity optimisation.")
            combiner = SingleDatasetIntensityCombiner(self, use_Imid)
        except DialsMergingStatisticsError as e:
            logger.info("Intensity combination failed with the error %s", e)
        else:
            intensity, variance = combiner.calculate_suitable_combined_intensities()
            # update data in reflection table
            isel = self.suitable_refl_for_scaling_sel.iselection()
            self._reflection_table["intensity"].set_selected(isel, intensity)
            self._reflection_table["variance"].set_selected(isel, variance)
            if self.global_Ih_table:
                self.global_Ih_table.update_data_in_blocks(
                    intensity, 0, column="intensity"
                )
                self.global_Ih_table.update_data_in_blocks(
                    variance, 0, column="variance"
                )
                self.global_Ih_table.update_weights(
                    self._experiment.scaling_model.error_model
                )
                self.global_Ih_table.calc_Ih()
            if self._free_Ih_table:
                self._free_Ih_table.update_data_in_blocks(
                    intensity, 0, column="intensity"
                )
                self._free_Ih_table.update_data_in_blocks(
                    variance, 0, column="variance"
                )
                self._free_Ih_table.update_weights(
                    self._experiment.scaling_model.error_model
                )
                self._free_Ih_table.calc_Ih()
            self.experiment.scaling_model.record_intensity_combination_Imid(
                combiner.max_key
            )

    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
        """
        Calculate scale factors for all suitable reflections.

        Use the current model to calculate scale factors for all suitable
        reflections, and set these in the reflection table. If caller=None,
        the global_Ih_table is updated. If calc_cov, an error estimate on the
        inverse scales is calculated.
        """
        self._reflection_table["inverse_scale_factor_variance"] = flex.double(
            self.reflection_table.size(), 0.0
        )
        n_blocks = self.params.scaling_options.nproc
        n_start = 0
        all_scales = flex.double([])
        all_invsfvars = flex.double([])
        n_param_tot = sum(c.n_params for c in self.components.values())
        for i in range(1, n_blocks + 1):  # do calc in blocks for speed/memory
            n_end = int(i * self.n_suitable_refl / n_blocks)
            block_isel = flex.size_t(range(n_start, n_end))
            n_start = n_end
            scales = flex.double(block_isel.size(), 1.0)
            scales_list = []
            derivs_list = []
            jacobian = sparse.matrix(block_isel.size(), n_param_tot)
            for component in self.components.values():
                component.update_reflection_data(block_selections=[block_isel])
                comp_scales, d = component.calculate_scales_and_derivatives(block_id=0)
                scales_list.append(comp_scales)
                if calc_cov:
                    derivs_list.append(d)
                scales *= comp_scales
            all_scales.extend(scales)
            if calc_cov and self.var_cov_matrix.non_zeroes > 0:
                n_cumulative_param = 0
                for j, component in enumerate(self.components):
                    d_block = derivs_list[j]
                    n_param = self.components[component].n_params
                    for k, component_2 in enumerate(self.components):
                        if component_2 != component:
                            d_block = row_multiply(d_block, scales_list[k])
                    jacobian.assign_block(d_block, 0, n_cumulative_param)
                    n_cumulative_param += n_param
                all_invsfvars.extend(
                    cpp_calc_sigmasq(jacobian.transpose(), self._var_cov_matrix)
                )
        scaled_isel = self.suitable_refl_for_scaling_sel.iselection()
        self.reflection_table["inverse_scale_factor"].set_selected(
            scaled_isel, all_scales
        )
        if calc_cov and self.var_cov_matrix.non_zeroes > 0:
            self.reflection_table["inverse_scale_factor_variance"].set_selected(
                scaled_isel, all_invsfvars
            )
        if caller is None:
            self.global_Ih_table.update_data_in_blocks(
                self.reflection_table["inverse_scale_factor"].select(
                    self.suitable_refl_for_scaling_sel
                ),
                dataset_id=0,
                column="inverse_scale_factor",
            )
            self.global_Ih_table.calc_Ih()
            if self._free_Ih_table:
                self._free_Ih_table.update_data_in_blocks(
                    self.reflection_table["inverse_scale_factor"].select(
                        self.suitable_refl_for_scaling_sel
                    ),
                    dataset_id=0,
                    column="inverse_scale_factor",
                )
                self._free_Ih_table.calc_Ih()
            logger.info(
                "Scale factors determined during minimisation have now been\n"
                "applied to all reflections for dataset %s.\n",
                self.reflection_table["id"][0],
            )

    def _select_reflections_for_scaling(self):
        """Select a subset of reflections to use in minimisation."""
        # For single dataset, auto means random
        if self.params.reflection_selection.method in (
            None,
            Auto,
            "auto",
            "random",
            "quasi_random",
        ):
            # do the random reflection selection.
            block = self.global_Ih_table.Ih_table_blocks[0]
            suitable_table = self.get_valid_reflections()
            presel = flumpy.to_numpy(
                calculate_scaling_subset_ranges(
                    suitable_table, self.params, print_summary=True
                )
            )
            preselection = presel[block.Ih_table["loc_indices"].to_numpy()]
            presel_block = block.select(preselection)
            # then select random groups
            n_h_over_1 = presel_block.calc_nh() > 1
            if n_h_over_1.any():
                presel_block = presel_block.select(n_h_over_1)
                min_groups = self.params.reflection_selection.random.min_groups
                avg_multi = np.mean(presel_block.group_multiplicities())
                min_refl = self.params.reflection_selection.random.min_reflections
                n_groups_in_table = presel_block.n_groups
                n_groups_to_sel = max(int(ceil(min_refl / avg_multi)), min_groups)

                logger.debug(
                    "Average multiplicity that reflection selection is sampling from: %s",
                    avg_multi,
                )
                if n_groups_to_sel < n_groups_in_table:
                    isel = flumpy.to_numpy(
                        flex.random_selection(n_groups_in_table, n_groups_to_sel)
                    )
                    loc_indices = (
                        presel_block.select_on_groups(isel)
                        .Ih_table["loc_indices"]
                        .to_numpy()
                    )
                else:
                    loc_indices = presel_block.Ih_table["loc_indices"].to_numpy()
                    n_groups_to_sel = n_groups_in_table
                self.scaling_selection = np.full(
                    self.n_suitable_refl, False, dtype=bool
                )
                # flex.bool(self.n_suitable_refl, False)
                self.scaling_selection[loc_indices] = True
                self.scaling_selection = flumpy.from_numpy(self.scaling_selection)
                logger.info(
                    """Randomly selected %s/%s groups (m>1) to use for scaling model
minimisation (%s reflections)""",
                    n_groups_to_sel,
                    n_groups_in_table,
                    loc_indices.size,
                )
            else:
                logger.warning(
                    """No groups left with multiplicity >1 after prefiltering,
attempting to use all reflections for minimisation."""
                )
                self.params.reflection_selection.method = "use_all"

        if self.params.reflection_selection.method == "intensity_ranges":
            overall_scaling_selection = calculate_scaling_subset_ranges_with_E2(
                self.reflection_table, self.params
            )
            self.scaling_selection = overall_scaling_selection.select(
                self.suitable_refl_for_scaling_sel
            )
            if self._free_Ih_table:
                self.scaling_selection.set_selected(self.free_set_selection, False)
        if self.params.reflection_selection.method == "use_all":
            if self._free_Ih_table:
                self.scaling_selection = ~self.free_set_selection
            else:
                self.scaling_selection = flex.bool(self.n_suitable_refl, True)
        self.scaling_subset_sel = copy.deepcopy(self.scaling_selection)
        self.scaling_selection &= ~self.outliers  # now apply outliers

    def _update_model_data(self):
        """Use the data in the Ih_table to update the model data."""
        assert self.Ih_table is not None
        block_selections = self.Ih_table.get_block_selections_for_dataset(dataset=0)
        block_selections = [flumpy.from_numpy(b) for b in block_selections]
        for component in self.components.values():
            component.update_reflection_data(block_selections=block_selections)

    def _create_Ih_table(self):
        """Create an Ih_table from the reflection table using the scaling selection."""
        self._Ih_table = IhTable(
            [self.get_reflections_for_model_minimisation()],
            self.experiment.crystal.get_space_group(),
            indices_lists=[self.scaling_selection.iselection()],
            nblocks=self.params.scaling_options.nproc,
            anomalous=self.params.anomalous,
        )
        if self._experiment.scaling_model.error_model:
            # update with the error model to add the correct weights
            self._Ih_table.update_weights(self._experiment.scaling_model.error_model)

    def _create_global_Ih_table(
        self, free_set_percentage=0, anomalous=False, remove_outliers=False
    ):
        sel_reflections = self.get_valid_reflections()
        if remove_outliers:
            sel_reflections = sel_reflections.select(~self.outliers)
        free_Ih_table = None
        global_Ih_table = IhTable(
            [sel_reflections],
            self.experiment.crystal.get_space_group(),
            nblocks=1,
            free_set_percentage=free_set_percentage,
            free_set_offset=self.params.scaling_options.free_set_offset,
            anomalous=anomalous,
        )
        if free_set_percentage:
            loc_indices = (
                global_Ih_table.blocked_data_list[-1].Ih_table["loc_indices"].to_numpy()
            )
            self.free_set_selection = flex.bool(self.n_suitable_refl, False)
            self.free_set_selection.set_selected(flumpy.from_numpy(loc_indices), True)
            global_Ih_table = IhTable(
                [sel_reflections.select(~self.free_set_selection)],
                self.experiment.crystal.get_space_group(),
                [(~self.free_set_selection).iselection()],
                nblocks=1,
                anomalous=anomalous,
            )
            free_Ih_table = IhTable(
                [sel_reflections.select(self.free_set_selection)],
                self.experiment.crystal.get_space_group(),
                [self.free_set_selection.iselection()],
                nblocks=1,
                anomalous=anomalous,
            )
        return global_Ih_table, free_Ih_table

    def _configure_model_and_datastructures(self, for_multi=False):
        """
        Store the relevant data in the scaling model components.

        This takes the columns from the 'suitable' part of the reflection table (
        which will include outliers). Then a global_Ih_table is created, which can
        be used for outlier rejection. When calculations are done with the model
        components, the correct reflections should first be selected out of the
        stored data.
        """
        sel_reflections = self.get_valid_reflections()
        self.experiment.scaling_model.configure_components(
            sel_reflections, self.experiment, self.params
        )
        free_set_percentage = 0
        if self.params.scaling_options.use_free_set and not for_multi:
            free_set_percentage = self.params.scaling_options.free_set_percentage
        self._global_Ih_table, self._free_Ih_table = self._create_global_Ih_table(
            free_set_percentage=free_set_percentage, anomalous=self.params.anomalous
        )
        rows = [[key, str(val.n_params)] for key, val in self.components.items()]
        logger.info("The following corrections will be applied to this dataset: \n")
        logger.info(tabulate(rows, ["correction", "n_parameters"]))

    @Subject.notify_event(event="performed_outlier_rejection")
    def round_of_outlier_rejection(self):
        """Perform a round of outlier rejection, set a new outliers array."""
        assert self.global_Ih_table is not None
        if self.params.scaling_options.outlier_rejection:
            outlier_indices = determine_outlier_index_arrays(
                self.global_Ih_table,
                self.params.scaling_options.outlier_rejection,
                self.params.scaling_options.outlier_zmax,
            )[0]
            outlier_indices = flumpy.from_numpy(outlier_indices)
            self.outliers = flex.bool(self.n_suitable_refl, False)
            self.outliers.set_selected(outlier_indices, True)
            if self.params.scaling_options.emax:
                Esq_outlier_indices = determine_Esq_outlier_index_arrays(
                    self.global_Ih_table,
                    self.experiment,
                    self.params.scaling_options.emax,
                )[0]
                self.outliers.set_selected(flumpy.from_numpy(Esq_outlier_indices), True)
            if self._free_Ih_table:
                free_outlier_indices = determine_outlier_index_arrays(
                    self._free_Ih_table,
                    self.params.scaling_options.outlier_rejection,
                    self.params.scaling_options.outlier_zmax,
                )[0]
                self.outliers.set_selected(
                    flumpy.from_numpy(free_outlier_indices), True
                )
                if self.params.scaling_options.emax:
                    free_Esq_outlier_indices = determine_Esq_outlier_index_arrays(
                        self._free_Ih_table,
                        self.experiment,
                        self.params.scaling_options.emax,
                    )[0]
                    self.outliers.set_selected(
                        flumpy.from_numpy(free_Esq_outlier_indices), True
                    )

    def make_ready_for_scaling(self, outlier=True):
        """
        Prepare the datastructures for a round of scaling.

        Update the scaling selection, create a new Ih_table and update the model
        data ready for minimisation.
        """
        if outlier:
            self.scaling_selection = self.scaling_subset_sel & ~self.outliers
        else:
            self.scaling_selection = copy.deepcopy(self.scaling_subset_sel)
        self._create_Ih_table()
        self._update_model_data()

    def clean_reflection_table(self):
        """Remove additional added columns that are not required for output."""
        self._initial_keys.append("inverse_scale_factor")
        self._initial_keys.append("inverse_scale_factor_variance")
        self._initial_keys.append("Ih_values")
        self._initial_keys.append("intensity.scale.value")
        self._initial_keys.append("intensity.scale.variance")
        self.reflection_table["intensity.scale.value"] = self.reflection_table[
            "intensity"
        ]
        self.reflection_table["intensity.scale.variance"] = self.reflection_table[
            "variance"
        ]
        if "Esq" in self.reflection_table:
            del self.reflection_table["Esq"]
        for key in self.reflection_table.keys():
            if key not in self._initial_keys:
                del self._reflection_table[key]
        bad = self._reflection_table.get_flags(
            self._reflection_table.flags.bad_for_scaling, all=False
        )
        self._reflection_table.set_flags(~bad, self.reflection_table.flags.scaled)


class MultiScalerBase(ScalerBase):
    """Base class for scalers handling multiple datasets."""

    def __init__(self, single_scalers):
        """Initialise from a list of single scalers."""
        super().__init__(single_scalers[0].params)
        self.single_scalers = single_scalers

    def remove_datasets(self, scalers, n_list):
        """
        Delete a scaler from the dataset.

        Code in this module does not necessarily have access to all references of
        experiments and reflections, so log the position in the list so that they
        can be deleted later. Scaling algorithm code should only depends on the
        scalers.
        """
        initial_number = len(scalers)
        for n in n_list[::-1]:
            self._removed_datasets.append(scalers[n].experiment.identifier)
            del scalers[n]
        # if 0 in n_list:
        #    self._experiment = scalers[0].experiments
        assert len(scalers) == initial_number - len(n_list)
        logger.info("Removed datasets: %s", n_list)
        if not self.active_scalers:
            raise ValueError("No datasets remain after filtering out bad datasets.")

    def get_free_set_reflections(self):
        """Get all reflections in the free set if it exists."""
        if self._free_Ih_table:
            refls = flex.reflection_table()
            for scaler in self.active_scalers:
                refls.extend(scaler.get_free_set_reflections())
            return refls
        return None

    def get_work_set_reflections(self):
        """Get all reflections in the free set if it exists."""
        if self._free_Ih_table:
            refls = flex.reflection_table()
            for scaler in self.active_scalers:
                refls.extend(scaler.get_work_set_reflections())
            return refls
        return None

    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
        """
        Calculate scale factors for all suitable reflections in the datasets.

        After the scale factors are updated, the global_Ih_table is updated also.
        """
        if calc_cov:
            logger.info("Calculating error estimates of inverse scale factors. \n")
        for i, scaler in enumerate(self.active_scalers):
            scaler.expand_scales_to_all_reflections(caller=self, calc_cov=calc_cov)
            # now update global Ih table
            self.global_Ih_table.update_data_in_blocks(
                scaler.reflection_table["inverse_scale_factor"].select(
                    scaler.suitable_refl_for_scaling_sel
                ),
                dataset_id=i,
                column="inverse_scale_factor",
            )
            if self._free_Ih_table:
                self._free_Ih_table.update_data_in_blocks(
                    scaler.reflection_table["inverse_scale_factor"].select(
                        scaler.suitable_refl_for_scaling_sel
                    ),
                    dataset_id=i,
                    column="inverse_scale_factor",
                )

        self.global_Ih_table.calc_Ih()
        if self._free_Ih_table:
            self._free_Ih_table.calc_Ih()
        logger.info(
            "Scale factors determined during minimisation have now been\n"
            "applied to all datasets.\n"
        )

    def update_for_minimisation(self, apm, block_id):
        """Update the scale factors and Ih for the next iteration of minimisation."""
        self._update_for_minimisation(apm, block_id, calc_Ih=True)

    def _update_for_minimisation(self, apm, block_id, calc_Ih=True):
        scales = flex.double([])
        derivs = []
        for apm_i in apm.apm_list:
            scales_i, derivs_i = RefinerCalculator.calculate_scales_and_derivatives(
                apm_i, block_id
            )
            scales.extend(scales_i)
            derivs.append(derivs_i)
        deriv_matrix = sparse.matrix(scales.size(), apm.n_active_params)
        start_row_no = 0
        for j, deriv in enumerate(derivs):
            deriv_matrix.assign_block(deriv, start_row_no, apm.apm_data[j]["start_idx"])
            start_row_no += deriv.n_rows
        self.Ih_table.set_inverse_scale_factors(flumpy.to_numpy(scales), block_id)
        self.Ih_table.set_derivatives(deriv_matrix, block_id)
        if calc_Ih:
            self.Ih_table.calc_Ih(block_id)
        # The parallelisation below would work if sparse matrices were
        # pickleable (I think!) - with more benefit for larger number of datasets.'''
        """def task_wrapper(block):
            s, d = RefinerCalculator.calculate_scales_and_derivatives(block)
            return s, d
        blocks = apm.apm_list
        task_results = easy_mp.parallel_map(func=task_wrapper, iterable=blocks,
            processes=n_datasets, method="multiprocessing",
            preserve_exception_message=True
        )
        scales_list, derivs_list = zip(*task_results)"""

    def _update_model_data(self):
        for i, scaler in enumerate(self.active_scalers):
            block_selections = self.Ih_table.get_block_selections_for_dataset(i)
            these_block_selections = [flumpy.from_numpy(b) for b in block_selections]
            for component in scaler.components.values():
                component.update_reflection_data(
                    block_selections=these_block_selections
                )

    def _create_global_Ih_table(self, anomalous=False, remove_outliers=False):
        if remove_outliers:
            tables = [
                s.get_valid_reflections().select(~s.outliers)
                for s in self.active_scalers
            ]
        else:
            tables = [s.get_valid_reflections() for s in self.active_scalers]
        free_set_percentage = 0.0
        free_Ih_table = None
        if self.params.scaling_options.use_free_set:
            free_set_percentage = self.params.scaling_options.free_set_percentage
        space_group = self.active_scalers[0].experiment.crystal.get_space_group()
        global_Ih_table = IhTable(
            tables,
            space_group,
            nblocks=1,
            additional_cols=["partiality"],
            free_set_percentage=free_set_percentage,
            free_set_offset=self.params.scaling_options.free_set_offset,
            anomalous=anomalous,
        )
        if free_set_percentage:
            # need to set free_set_selection in individual scalers
            tables = []
            free_tables = []
            indices_list = []
            free_indices_list = []
            for i, scaler in enumerate(self.active_scalers):
                sel = (
                    global_Ih_table.Ih_table_blocks[-1]
                    .Ih_table["dataset_id"]
                    .to_numpy()
                    == i
                )
                indiv_Ih_block = global_Ih_table.Ih_table_blocks[-1].select(sel)
                loc_indices = indiv_Ih_block.Ih_table["loc_indices"].to_numpy()
                scaler.free_set_selection = flex.bool(scaler.n_suitable_refl, False)
                scaler.free_set_selection.set_selected(
                    flumpy.from_numpy(loc_indices), True
                )
                tables.append(scaler.get_work_set_reflections())
                free_tables.append(scaler.get_free_set_reflections())
                free_indices_list.append(scaler.free_set_selection.iselection())
                indices_list.append((~scaler.free_set_selection).iselection())
            global_Ih_table = IhTable(
                tables,
                space_group,
                indices_list,
                nblocks=1,
                anomalous=anomalous,
            )
            free_Ih_table = IhTable(
                free_tables,
                space_group,
                free_indices_list,
                nblocks=1,
                anomalous=anomalous,
            )
        return global_Ih_table, free_Ih_table

    def _create_Ih_table(self):
        """Create a new Ih table from the reflection tables."""
        tables = [
            s.get_reflections_for_model_minimisation() for s in self.active_scalers
        ]
        indices_lists = [s.scaling_selection.iselection() for s in self.active_scalers]
        self._Ih_table = IhTable(
            tables,
            self.active_scalers[0].experiment.crystal.get_space_group(),
            indices_lists=indices_lists,
            nblocks=self.params.scaling_options.nproc,
            anomalous=self.params.anomalous,
        )
        for i, scaler in enumerate(self.active_scalers):
            error_model = scaler._experiment.scaling_model.error_model
            if error_model:
                self._Ih_table.update_weights(error_model, dataset_id=i)

    def make_ready_for_scaling(self, outlier=True):
        """
        Prepare the datastructures for a round of scaling.

        Update the scaling selection, create a new Ih_table and update the model
        data ready for minimisation. Also check to see if any datasets should be
        removed.
        """
        datasets_to_remove = []
        for i, scaler in enumerate(self.active_scalers):
            if outlier:
                scaler.scaling_selection = scaler.scaling_subset_sel & ~scaler.outliers
            else:
                scaler.scaling_selection = copy.deepcopy(scaler.scaling_subset_sel)
            if scaler.scaling_selection.count(True) == 0:
                datasets_to_remove.append(i)
        if datasets_to_remove:
            self.remove_datasets(self.active_scalers, datasets_to_remove)
            (
                self._global_Ih_table,
                self._free_Ih_table,
            ) = self._create_global_Ih_table(self.params.anomalous)
        self._create_Ih_table()
        self._update_model_data()

    @Subject.notify_event(event="performed_outlier_rejection")
    def round_of_outlier_rejection(self):
        self._round_of_outlier_rejection(target=None)

    def _round_of_outlier_rejection(self, target=None):
        """
        Perform a round of outlier rejection across all datasets.

        After identifying outliers, set the outliers property in individual scalers.
        """
        assert self.active_scalers is not None
        if not self.global_Ih_table:
            self._global_Ih_table, self._free_Ih_table = self._create_global_Ih_table(
                anomalous=self.params.anomalous
            )
        if self.params.scaling_options.outlier_rejection:
            if target:
                method = "target"
            else:
                method = self.params.scaling_options.outlier_rejection
            outlier_index_arrays = determine_outlier_index_arrays(
                self.global_Ih_table,
                method,
                self.params.scaling_options.outlier_zmax,
                target=target,
            )
            for outlier_indices, scaler in zip(
                outlier_index_arrays, self.active_scalers
            ):
                scaler.outliers = flex.bool(scaler.n_suitable_refl, False)
                scaler.outliers.set_selected(flumpy.from_numpy(outlier_indices), True)
            if self.params.scaling_options.emax:
                Esq_outlier_indices = determine_Esq_outlier_index_arrays(
                    self.global_Ih_table,
                    self.active_scalers[0].experiment,
                    self.params.scaling_options.emax,
                )
                for Esq_outliers, scaler in zip(
                    Esq_outlier_indices, self.active_scalers
                ):
                    scaler.outliers.set_selected(flumpy.from_numpy(Esq_outliers), True)
            if self._free_Ih_table:
                free_outlier_index_arrays = determine_outlier_index_arrays(
                    self._free_Ih_table,
                    method,
                    self.params.scaling_options.outlier_zmax,
                    target=target,
                )
                for outlier_indices, scaler in zip(
                    free_outlier_index_arrays,
                    self.active_scalers,
                ):
                    scaler.outliers.set_selected(
                        flumpy.from_numpy(outlier_indices), True
                    )
                if self.params.scaling_options.emax:
                    free_Esq_outlier_indices = determine_Esq_outlier_index_arrays(
                        self._free_Ih_table,
                        self.active_scalers[0].experiment,
                        self.params.scaling_options.emax,
                    )
                    for Esq_outlier_indices, scaler in zip(
                        free_Esq_outlier_indices,
                        self.active_scalers,
                    ):
                        scaler.outliers.set_selected(
                            flumpy.from_numpy(Esq_outlier_indices), True
                        )

        logger.debug("Finished outlier rejection.")
        log_memory_usage()

    def _select_reflections_for_scaling(self):
        # For multi dataset, auto means quasi random
        if self.params.reflection_selection.method in (
            None,
            Auto,
            "auto",
            "quasi_random",
        ):
            random_phil = self.params.reflection_selection.random
            block = self.global_Ih_table.Ih_table_blocks[0]
            m = block.group_multiplicities()
            mover1 = m[m > 1]
            total_refl_available = np.sum(mover1)
            avg_multi_overall = np.mean(mover1)
            n_datasets = len(self.active_scalers)
            total_target = max(
                avg_multi_overall * random_phil.min_groups, random_phil.min_reflections
            )

            # strategy - want ~1000 groups worth with quasi-random.
            target_for_qr = min(1000, random_phil.min_groups) * avg_multi_overall
            # ^ select this many for qr, then some from individual, then some
            # from overall.
            indices, dataset_ids = select_connected_reflections_across_datasets(
                self.global_Ih_table,
                self.active_scalers[0].experiment,
                random_phil.multi_dataset.Isigma_cutoff,
                min_total=target_for_qr,
            )

            n_connected_by_dataset = []
            total_qr_sel = indices.size()
            for i, scaler in enumerate(self.active_scalers):
                scaler.scaling_selection = flex.bool(scaler.n_suitable_refl, False)
                sel = dataset_ids == i
                indices_for_dataset = indices.select(sel)
                scaler.scaling_selection.set_selected(indices_for_dataset, True)
                n_connected_by_dataset.append(indices_for_dataset.size())

            min_cross_dataset = (total_target - total_qr_sel) / (2.0 * n_datasets)

            # first select individual dataset reflections, as we can then
            # calculate how many this procedure selects.
            total_individual_selection = 0
            total_indiv_dataset = flex.double()
            total_used_so_far = 0

            global_block = _select_groups_on_Isigma_cutoff(
                self.global_Ih_table.Ih_table_blocks[0],
                random_phil.multi_dataset.Isigma_cutoff,
            )
            for i, scaler in enumerate(self.active_scalers):
                # first set cross dataset connected reflections identified.
                # scaler.scaling_selection = flex.bool(scaler.n_suitable_refl, False)
                # now select randomly within dataset, first do preselection.
                if min_cross_dataset > 0:
                    sel = global_block.Ih_table["dataset_id"].to_numpy() == i
                    indiv_Ih_block = global_block.select(sel)
                    suitable_table = scaler.reflection_table.select(
                        scaler.suitable_refl_for_scaling_sel
                    )
                    presel = flumpy.to_numpy(
                        calculate_scaling_subset_ranges(suitable_table, self.params)
                    )
                    preselection = presel[
                        indiv_Ih_block.Ih_table["loc_indices"].to_numpy()
                    ]
                    block = indiv_Ih_block.select(preselection)
                    # have a block of individual dataset to chose from.

                    # now select multi > 1
                    n_h_over_1 = block.calc_nh() > 1
                    block = block.select(n_h_over_1)
                    if block.size:
                        n_groups = int(
                            ceil(
                                min_cross_dataset
                                / np.mean(block.group_multiplicities())
                            )
                        )
                        n_groups_in_table = block.n_groups
                        if n_groups < n_groups_in_table:
                            isel = flumpy.to_numpy(
                                flex.random_selection(n_groups_in_table, n_groups)
                            )
                            loc_indices = (
                                block.select_on_groups(isel)
                                .Ih_table["loc_indices"]
                                .to_numpy()
                            )
                            scaler.scaling_selection.set_selected(
                                flumpy.from_numpy(loc_indices), True
                            )
                        else:
                            loc_indices = block.Ih_table["loc_indices"].to_numpy()
                            scaler.scaling_selection.set_selected(
                                flumpy.from_numpy(loc_indices), True
                            )
                        n_refl = loc_indices.size
                    else:  # no refl with multiplicity > 1.
                        n_refl = 0
                else:
                    n_refl = 0
                total_individual_selection += n_refl
                total_indiv_dataset.append(n_refl)
                total_used_so_far += scaler.scaling_selection.count(True)

            # now do random selection across datasets
            target_left = total_target - total_used_so_far
            pc_used = total_used_so_far / total_refl_available
            if target_left > 0 and pc_used < 1:
                min_refl_needed = target_left * 1.05 / (1.0 - pc_used)
                # now do random sel on all
                n_groups_to_sel = int(min_refl_needed / avg_multi_overall)

                n_groups_in_table = global_block.n_groups
                if n_groups_to_sel < n_groups_in_table:
                    isel = flumpy.to_numpy(
                        flex.random_selection(n_groups_in_table, n_groups_to_sel)
                    )
                    sel_block = global_block.select_on_groups(isel)
                else:  # just use all
                    sel_block = global_block
            else:
                sel_block = global_block.select_on_groups(
                    np.array([], dtype=np.uint64).reshape((0,))
                )

            header = [
                "Dataset id",
                "reflections \nconnected to \nother datasets",
                "randomly selected \nreflections \nwithin dataset",
                "randomly selected \nreflections \nacross datasets",
                "combined number \nof reflections",
            ]
            rows = []
            total_overall = 0

            for i, scaler in enumerate(self.active_scalers):
                sel = sel_block.Ih_table["dataset_id"].to_numpy() == i
                indiv_block = sel_block.select(sel)
                loc_indices = indiv_block.Ih_table["loc_indices"].to_numpy()
                scaler.scaling_selection.set_selected(
                    flumpy.from_numpy(loc_indices), True
                )
                n = scaler.scaling_selection.count(True)
                rows.append(
                    [
                        scaler.reflection_table.experiment_identifiers().keys()[0],
                        str(n_connected_by_dataset[i]),
                        str(total_indiv_dataset[i]),
                        str(loc_indices.size),
                        str(n),
                    ]
                )
                total_overall += n
                scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                scaler.scaling_selection &= ~scaler.outliers

            rows.append(
                [
                    "total",
                    str(total_qr_sel),
                    str(total_individual_selection),
                    str(sel_block.size),
                    str(total_overall),
                ]
            )
            logger.info(
                "Summary of reflections chosen for minimisation from each dataset (%s total):",
                total_overall,
            )
            logger.info(tabulate(rows, header))

        elif self.params.reflection_selection.method == "intensity_ranges":
            datasets_to_remove = []
            for i, scaler in enumerate(self.active_scalers):
                try:
                    overall_scaling_selection = calculate_scaling_subset_ranges_with_E2(
                        scaler.reflection_table, scaler.params
                    )
                except BadDatasetForScalingException:
                    datasets_to_remove.append(i)
                else:
                    scaler.scaling_selection = overall_scaling_selection.select(
                        scaler.suitable_refl_for_scaling_sel
                    )
                    if self._free_Ih_table:
                        scaler.scaling_selection.set_selected(
                            scaler.free_set_selection, False
                        )
                    scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                    scaler.scaling_selection &= ~scaler.outliers
            if datasets_to_remove:
                self.remove_datasets(self.active_scalers, datasets_to_remove)
                # Reinitialise the global datastructures for the reduced dataset.
                (
                    self._global_Ih_table,
                    self._free_Ih_table,
                ) = self._create_global_Ih_table(self.params.anomalous)
        elif self.params.reflection_selection.method == "use_all":
            block = self.global_Ih_table.Ih_table_blocks[0]
            n_mult = np.count_nonzero(block.calc_nh() > 1)
            logger.info(
                "Using all reflections (%s) for minimisation (%s with m>1)",
                block.size,
                n_mult,
            )
            for scaler in self.active_scalers:
                if self._free_Ih_table:
                    scaler.scaling_selection = ~scaler.free_set_selection
                else:
                    scaler.scaling_selection = flex.bool(scaler.n_suitable_refl, True)
                scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                scaler.scaling_selection &= ~scaler.outliers
        elif self.params.reflection_selection.method == "random":
            # random means a random subset of groups,
            random_phil = self.params.reflection_selection.random
            block = self.global_Ih_table.Ih_table_blocks[0]
            sel_block = block.select(block.calc_nh() > 1)

            block = _select_groups_on_Isigma_cutoff(
                block, random_phil.multi_dataset.Isigma_cutoff
            )
            block = block.select(block.calc_nh() > 1)
            if not block.size:
                raise SystemExit(
                    "No groups left with multiplicity >1, scaling not possible."
                )
            avg_multi = np.mean(block.group_multiplicities())
            n_groups_to_sel = max(
                random_phil.min_groups,
                int(ceil(random_phil.min_reflections / avg_multi)),
            )

            n_groups_in_table = block.n_groups
            if n_groups_to_sel < n_groups_in_table:
                isel = flex.random_selection(n_groups_in_table, n_groups_to_sel)
                sel_block = block.select_on_groups(flumpy.to_numpy(isel))
            else:  # just use all
                sel_block = block
            logger.info(
                "Selected %s/%s groups (m>1) to use for minimisation (%s reflections)",
                min(n_groups_to_sel, n_groups_in_table),
                n_groups_in_table,
                sel_block.size,
            )
            for i, scaler in enumerate(self.active_scalers):
                sel = sel_block.Ih_table["dataset_id"].to_numpy() == i
                indiv_block = sel_block.select(sel)
                loc_indices = indiv_block.Ih_table["loc_indices"].to_numpy()
                scaler.scaling_selection = flex.bool(scaler.n_suitable_refl, False)
                scaler.scaling_selection.set_selected(
                    flumpy.from_numpy(loc_indices), True
                )
                scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                scaler.scaling_selection &= ~scaler.outliers
        else:
            raise ValueError("Invalid choice for 'reflection_selection.method'.")

    def determine_shared_model_components(self):
        shared_components = set()
        for scaler in self.active_scalers:
            shared_components.add(
                scaler.experiment.scaling_model.get_shared_components()
            )
        if len(shared_components) == 1:
            return list(shared_components)[0]
        return None

    @Subject.notify_event(event="performed_error_analysis")
    def perform_error_optimisation(self, update_Ih=True):
        """Perform an optimisation of the sigma values."""
        if self.params.weighting.error_model.grouping == "combined":
            minimisation_groups = [[i for i, _ in enumerate(self.active_scalers)]]
        elif self.params.weighting.error_model.grouping == "individual":
            minimisation_groups = [[i] for i, _ in enumerate(self.active_scalers)]
        else:
            groups = self.params.weighting.error_model.error_model_group
            if not groups:
                logger.info(
                    """No error model groups defined, defaulting to combined error model optimisation"""
                )
                minimisation_groups = [[i for i, _ in enumerate(self.active_scalers)]]
            else:
                all_datasets = [i for i, _ in enumerate(self.active_scalers)]
                # groups are defined in terms of sweeps (1,2,3,...), but here
                # need to convert to dataset number (0, 1, 2,...)
                explicitly_grouped = [i - 1 for j in groups for i in j]
                if -1 in explicitly_grouped:  # sweeps provided indexed from 0
                    explicitly_grouped = [i for j in groups for i in j]
                    minimisation_groups = [list(g) for g in groups]
                else:
                    minimisation_groups = [[i - 1 for i in g] for g in groups]
                others = set(all_datasets).difference(set(explicitly_grouped))
                if others:
                    minimisation_groups += [list(others)]

        for g in minimisation_groups:
            scalers = [self.active_scalers[i] for i in g]
            error_model = scalers[0]._experiment.scaling_model.error_model
            if not error_model.params.minimisation:
                continue
            tables = [s.get_valid_reflections().select(~s.outliers) for s in scalers]
            space_group = scalers[0].experiment.crystal.get_space_group()
            Ih_table = IhTable(tables, space_group, anomalous=True)
            if len(minimisation_groups) == 1:
                logger.info("Determining a combined error model for all datasets")
            else:
                logger.info(
                    f"Error model determination for sweep(s) {','.join(str(i+1) for i in g)}"
                )
            try:
                model = run_error_model_refinement(
                    scalers[0]._experiment.scaling_model.error_model, Ih_table
                )
            except (ValueError, RuntimeError) as e:
                logger.info(e)
                logger.debug(e, exc_info=True)
            else:
                for s in scalers:
                    s._experiment.scaling_model.set_error_model(model)
        if update_Ih:
            if len(minimisation_groups) == 1:
                # one overall error model was refined
                error_model = self.active_scalers[
                    0
                ]._experiment.scaling_model.error_model
                if not error_model.params.minimisation:
                    return  # no need to update
                self.global_Ih_table.update_weights(error_model)
                if self._free_Ih_table:
                    self._free_Ih_table.update_weights(error_model)
            else:
                for i, scaler in enumerate(self.active_scalers):
                    error_model = scaler._experiment.scaling_model.error_model
                    if not error_model.params.minimisation:
                        continue  # no need to update for this subset
                    self.global_Ih_table.update_weights(error_model, dataset_id=i)
                    if self._free_Ih_table:
                        self._free_Ih_table.update_weights(error_model, dataset_id=i)


class MultiScaler(MultiScalerBase):
    """Scaler for multiple datasets where all datasets are being minimised."""

    id_ = "multi"

    def __init__(self, single_scalers):
        """
        Initialise a multiscaler from a list of single scalers.

        Create a global_Ih_table, an Ih_table to use for minimisation and update
        the data in the model components.
        """
        logger.info("Configuring a MultiScaler to handle the individual Scalers. \n")
        super().__init__(single_scalers)
        logger.info("Determining symmetry equivalent reflections across datasets.\n")
        self._active_scalers = self.single_scalers
        self._global_Ih_table, self._free_Ih_table = self._create_global_Ih_table(
            self.params.anomalous
        )
        # now select reflections from across the datasets
        self._select_reflections_for_scaling()
        self._create_Ih_table()
        # now add data to scale components from datasets
        self._update_model_data()
        logger.info("Completed configuration of MultiScaler. \n\n" + "=" * 80 + "\n")
        log_memory_usage()

    def fix_initial_parameter(self):
        for scaler in self.active_scalers:
            fixed = scaler.experiment.scaling_model.fix_initial_parameter(self.params)
            if fixed:
                return fixed

    def combine_intensities(self):
        """Combine reflection intensities, either jointly or separately."""
        if self.params.reflection_selection.combine.joint_analysis:
            logger.info(
                "Performing multi-dataset profile/summation intensity optimisation."
            )
            try:
                combiner = MultiDatasetIntensityCombiner(self)
            except DialsMergingStatisticsError as e:
                logger.info("Intensity combination failed with the error %s", e)
            else:
                for i, scaler in enumerate(self.active_scalers):
                    (
                        intensity,
                        variance,
                    ) = combiner.calculate_suitable_combined_intensities(i)
                    isel = scaler.suitable_refl_for_scaling_sel.iselection()
                    scaler.reflection_table["intensity"].set_selected(isel, intensity)
                    scaler.reflection_table["variance"].set_selected(isel, variance)
                    self.global_Ih_table.update_data_in_blocks(
                        intensity, i, column="intensity"
                    )
                    self.global_Ih_table.update_data_in_blocks(
                        variance, i, column="variance"
                    )
                    self.global_Ih_table.update_weights(
                        error_model=scaler._experiment.scaling_model.error_model,
                        dataset_id=i,
                    )
                    scaler.experiment.scaling_model.record_intensity_combination_Imid(
                        combiner.max_key
                    )
                    if self._free_Ih_table:
                        self._free_Ih_table.update_data_in_blocks(
                            intensity, i, column="intensity"
                        )
                        self._free_Ih_table.update_data_in_blocks(
                            variance, i, column="variance"
                        )
                        self._free_Ih_table.update_weights(
                            error_model=scaler._experiment.scaling_model.error_model,
                            dataset_id=i,
                        )
                self.global_Ih_table.calc_Ih()
                if self._free_Ih_table:
                    self._free_Ih_table.calc_Ih()
        else:
            for i, scaler in enumerate(self.single_scalers):
                scaler.combine_intensities()
                intensity = scaler.reflection_table["intensity"].select(
                    scaler.suitable_refl_for_scaling_sel
                )
                variance = scaler.reflection_table["variance"].select(
                    scaler.suitable_refl_for_scaling_sel
                )
                self.global_Ih_table.update_data_in_blocks(
                    intensity, i, column="intensity"
                )
                self.global_Ih_table.update_data_in_blocks(
                    variance, i, column="variance"
                )
                self.global_Ih_table.update_weights(
                    error_model=scaler._experiment.scaling_model.error_model,
                    dataset_id=i,
                )
                if self._free_Ih_table:
                    self._free_Ih_table.update_data_in_blocks(
                        intensity, i, column="intensity"
                    )
                    self._free_Ih_table.update_data_in_blocks(
                        variance, i, column="variance"
                    )
                    self._free_Ih_table.update_weights(
                        error_model=scaler._experiment.scaling_model.error_model,
                        dataset_id=i,
                    )
            self.global_Ih_table.calc_Ih()
            if self._free_Ih_table:
                self._free_Ih_table.calc_Ih()


class TargetScaler(MultiScalerBase):
    """A target scaler for scaling datasets against already scaled data."""

    id_ = "target"

    def __init__(self, scaled_scalers, unscaled_scalers):
        """
        Initialise a multiscaler from a list of single and unscaled scalers.

        First, set the active scalers (the unscaled scalers) and use these to
        create a global_Ih_table. Then, use the scaled_scalers to create a
        target_Ih_table. Create an Ih_table to use for minimisation and use
        the target_Ih_table to set the Ih_values. Finally, update the data in
        the model components.
        """
        logger.info("\nInitialising a TargetScaler instance. \n")
        super().__init__(scaled_scalers)
        logger.info("Determining symmetry equivalent reflections across datasets.\n")
        self.unscaled_scalers = unscaled_scalers
        self._active_scalers = unscaled_scalers
        self._target_active_scalers = unscaled_scalers
        self._global_Ih_table, self._free_Ih_table = self._create_global_Ih_table(
            self.params.anomalous
        )
        self._select_reflections_for_scaling()
        tables = [
            s.reflection_table.select(s.suitable_refl_for_scaling_sel).select(
                s.scaling_selection
            )
            for s in self.single_scalers
        ]
        self._target_Ih_table = IhTable(
            tables,
            self.active_scalers[0].experiment.crystal.get_space_group(),
            nblocks=1,
            anomalous=self.params.anomalous,
        )  # Keep in one table for matching below
        self._create_Ih_table()
        self._update_model_data()
        logger.info("Completed initialisation of TargetScaler. \n\n" + "=" * 80 + "\n")
        log_memory_usage()

    @property
    def target_Ih_table(self):
        """An Ih_table containing data for the target."""
        return self._target_Ih_table

    def _create_Ih_table(self):
        super()._create_Ih_table()
        for block in self._Ih_table.blocked_data_list:
            # this step reduces the number of reflections in each block
            block.match_Ih_values_to_target(self._target_Ih_table)
        self.Ih_table.generate_block_selections()

    def round_of_outlier_rejection(self):
        """Perform a round of targeted outlier rejection."""
        self._round_of_outlier_rejection(target=self._target_Ih_table)

    def update_for_minimisation(self, apm, block_id):
        """Calculate the new parameters but don't calculate a new Ih."""
        self._update_for_minimisation(apm, block_id, calc_Ih=False)

    @Subject.notify_event(event="performed_scaling")
    def perform_scaling(self, engine=None, max_iterations=None, tolerance=None):
        """Minimise the scaling model, using a fixed-Ih target."""
        for scaler in self._target_active_scalers:
            self._active_scalers = [scaler]
            self._create_Ih_table()
            self._update_model_data()
            self._perform_scaling(
                target_type=ScalingTargetFixedIH,
                engine=engine,
                max_iterations=max_iterations,
                tolerance=tolerance,
            )
        self._active_scalers = self._target_active_scalers


class NullScaler(ScalerBase):
    """A singlescaler to allow targeted scaling against calculated intensities."""

    id_ = "null"

    def __init__(self, params, experiment, reflection):
        """Set the required properties to use as a scaler for targeted scaling."""
        super().__init__(params)
        self._experiment = experiment
        self._reflection_table = reflection
        self.n_suitable_refl = self._reflection_table.size()
        self._reflection_table["inverse_scale_factor"] = flex.double(
            self.n_suitable_refl, 1.0
        )
        if "variance" not in self._reflection_table:
            self._reflection_table["variance"] = flex.double(self.n_suitable_refl, 1.0)
        self._reflection_table.set_flags(
            flex.bool(self.n_suitable_refl, False),
            self._reflection_table.flags.excluded_for_scaling,
        )
        self.suitable_refl_for_scaling_sel = flex.bool(self.n_suitable_refl, True)
        self.outliers = flex.bool(self.n_suitable_refl, False)
        self.scaling_selection = flex.bool(self.n_suitable_refl, True)
        logger.info("Target dataset contains %s reflections", self.n_suitable_refl)
        logger.info(
            "Completed preprocessing and initialisation for this dataset."
            "\n\n" + "=" * 80 + "\n"
        )

    @property
    def experiment(self):
        """Return the experiment object for the dataset"""
        return self._experiment

    @property
    def reflection_table(self):
        """Return the reflection_table object for the dataset"""
        return self._reflection_table

    @property
    def components(self):
        """Shortcut to scaling model components."""
        return self.experiment.scaling_model.components

    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
        """Fill in abstract method, do nothing."""

    def update_for_minimisation(self, apm, block_id=0):
        """Fill in abstract method, do nothing."""


def calc_sf_variances(components, var_cov):
    """Calculate the variances of the inverse scales."""
    # note - can we do this calculation blockwise as well - takes quite a bit of memory?
    n_param = 0
    for component in components:
        n_param += components[component].n_params
        n_refl = sum(components[component].n_refl)  # should all be same
    jacobian = sparse.matrix(n_refl, n_param)
    n_cumulative_param = 0
    scales_list = []
    derivs_list = []
    for component in components:
        s, d = components[component].calculate_scales_and_derivatives(block_id=0)
        scales_list.append(s)
        derivs_list.append(d)
    for i, component in enumerate(components):
        d_block = derivs_list[i]
        n_param = components[component].n_params
        for j, component_2 in enumerate(components):
            if component_2 != component:
                d_block = row_multiply(d_block, scales_list[j])
        jacobian.assign_block(d_block, 0, n_cumulative_param)
        n_cumulative_param += n_param
    return cpp_calc_sigmasq(jacobian.transpose(), var_cov)
