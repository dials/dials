"""
This module defines classes which implement the stages of the scaling algorithm.

These 'scalers' act to initialise and connect various parts of the scaling
algorithm and datastructures such as the Ih_table, basis_function etc, and
present a united interface to the main scaling algorithm for single, multi
and targeted scaling.

The SingleScaler is defined, for scaling of a single dataset, a MultiScaler is
defined for scaling multiple datasets simultaneously and a TargetScaler is
defined for targeted scaling.
"""
from __future__ import absolute_import, division, print_function

import abc
import copy
import logging
import time
from collections import OrderedDict

import six
from cctbx import crystal, sgtbx
from dials_scaling_ext import row_multiply
from dials_scaling_ext import calc_sigmasq as cpp_calc_sigmasq
from dials.array_family import flex
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.outlier_rejection import determine_outlier_index_arrays
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.target_function import ScalingTarget, ScalingTargetFixedIH
from dials.algorithms.scaling.scaling_refiner import (
    scaling_refinery,
    error_model_refinery,
)
from dials.algorithms.scaling.error_model.error_model import get_error_model
from dials.algorithms.scaling.error_model.error_model_target import ErrorModelTarget
from dials.algorithms.scaling.parameter_handler import create_apm_factory
from dials.algorithms.scaling.scaling_utilities import (
    log_memory_usage,
    DialsMergingStatisticsError,
)
from dials.algorithms.scaling.combine_intensities import (
    SingleDatasetIntensityCombiner,
    MultiDatasetIntensityCombiner,
)
from dials.algorithms.scaling.reflection_selection import (
    calculate_scaling_subset_connected,
    calculate_scaling_subset_ranges_with_E2,
    calculate_scaling_subset_ranges,
    select_connected_reflections_across_datasets,
)
from dials.util.observer import Subject
from libtbx.table_utils import simple_table
from scitbx import sparse

logger = logging.getLogger("dials")


class ScalerBase(Subject):
    """
    Abstract base class for all scalers (single and multiple).
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        """Define the properties of a scaler."""
        super(ScalerBase, self).__init__(
            events=[
                "performed_scaling",
                "performed_error_analysis",
                "performed_outlier_rejection",
            ]
        )
        self._experiment = None
        self._space_group = None
        self._params = None
        self._reflection_table = []
        self._Ih_table = None
        self._global_Ih_table = None
        self._initial_keys = []
        self._basis_function = basis_function()
        self._final_rmsds = []
        self._removed_datasets = []
        self.error_model = None

    @property
    def removed_datasets(self):
        """The list of removed datasets."""
        return self._removed_datasets

    @property
    def final_rmsds(self):
        """Holder for final R-factors from last minimisation."""
        return self._final_rmsds

    @final_rmsds.setter
    def final_rmsds(self, new_rmsds):
        self._final_rmsds = new_rmsds

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
    def experiment(self):
        """The experiment object associated with the dataset."""
        return self._experiment

    @property
    def space_group(self):
        """The space_group associated with the dataset."""
        return self._space_group

    @space_group.setter
    def space_group(self, new_sg):
        if isinstance(new_sg, str):
            crystal_symmetry = crystal.symmetry(space_group_symbol=new_sg)
            self._space_group = crystal_symmetry.space_group()
        elif isinstance(new_sg, sgtbx.space_group):
            self._space_group = new_sg
        else:
            raise AssertionError(
                """Space group not recognised as a space group symbol
        or cctbx.sgtbx.space group object."""
            )
        if self.experiment:
            self.experiment.crystal.set_space_group(self._space_group)

    @property
    def reflection_table(self):
        """The reflection table of the datatset."""
        return self._reflection_table

    @reflection_table.setter
    def reflection_table(self, new_table):
        """Set the reflection table of the datatset."""
        self._reflection_table = new_table

    @property
    def params(self):
        """The params phil scope."""
        return self._params

    @property
    def initial_keys(self):
        """A list of initial reflection table keys."""
        return self._initial_keys

    @abc.abstractmethod
    def update_for_minimisation(self, apm, block_id):
        """Update the scale factors and Ih for the next minimisation iteration."""

    @abc.abstractmethod
    def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
        """Expand scales from a subset to all reflections."""

    @Subject.notify_event(event="performed_scaling")
    def perform_scaling(
        self, target_type=ScalingTarget, engine=None, max_iterations=None
    ):
        """Minimise the scaling model."""
        apm_factory = create_apm_factory(self)
        for _ in range(apm_factory.n_cycles):
            apm = apm_factory.make_next_apm()
            if not engine:
                engine = self.params.scaling_refinery.engine
            if not max_iterations:
                max_iterations = self.params.scaling_refinery.max_iterations
            st = time.time()
            refinery = scaling_refinery(
                engine=engine,
                scaler=self,
                target=target_type(),
                prediction_parameterisation=apm,
                max_iterations=max_iterations,
            )
            try:
                refinery.run()
            except Exception as e:
                logger.error(e, exc_info=True)
            ft = time.time()
            logger.info("Time taken for refinement %s", (ft - st))
            refinery.return_scaler()
            logger.info("\n" + "=" * 80 + "\n")

    @Subject.notify_event(event="performed_error_analysis")
    def perform_error_optimisation(self, update_Ih=True):
        """Perform an optimisation of the sigma values."""
        Ih_table = self.global_Ih_table
        Ih_table.reset_error_model()
        Ih_table.calc_Ih()
        error_model = get_error_model(self.params.weighting.error_model.error_model)
        try:
            refinery = error_model_refinery(
                engine="SimpleLBFGS",
                target=ErrorModelTarget(
                    error_model(
                        Ih_table.blocked_data_list[0],
                        self.params.weighting.error_model.n_bins,
                        self.params.weighting.error_model.min_Ih,
                        self.params.reflection_selection.min_partiality,
                    )
                ),
                max_iterations=100,
            )
            refinery.run()
        except (RuntimeError, ValueError) as e:
            logger.error(e, exc_info=True)
        else:
            error_model = refinery.return_error_model()
            logger.info(error_model)
            error_model.minimisation_summary()
            self.update_error_model(error_model, update_Ih=update_Ih)
        return error_model

    def clear_memory_from_derivs(self, block_id):
        """Remove derivatives from Ih_table if no longer needed."""
        del self.Ih_table.blocked_data_list[block_id].derivatives

    def clear_Ih_table(self):
        """Delete the data from the current Ih_table."""
        self._Ih_table = []


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
        super(SingleScaler, self).__init__()
        self._experiment = experiment
        self._params = params
        self.active_scalers = [self]
        self._space_group = self.experiment.crystal.get_space_group()
        n_model_params = sum(val.n_params for val in self.components.values())
        self._var_cov = sparse.matrix(n_model_params, n_model_params)
        self._initial_keys = list(reflection_table.keys())
        self._reflection_table = reflection_table
        self._Ih_table = None  # stores data for reflections used for minimisation
        self.suitable_refl_for_scaling_sel = self.get_suitable_for_scaling_sel(
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
            None
        )  # A selection of len n_suitable_refl of scaling subset selection
        self.scaling_selection = None  # As above, but with outliers deselected also
        self._configure_model_and_datastructures()
        if "Imid" in self.experiment.scaling_model.configdict:
            self.combine_intensities(self.experiment.scaling_model.configdict["Imid"])
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

    @staticmethod
    def get_suitable_for_scaling_sel(reflections):
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
    def consecutive_refinement_order(self):
        """Link to consecutive refinement order for parameter manager."""
        return self.experiment.scaling_model.consecutive_refinement_order

    @property
    def var_cov_matrix(self):
        """The variance covariance matrix for the parameters."""
        return self._var_cov

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
            self._var_cov.assign_block(
                var_cov_list.matrix_copy_block(
                    0, 0, apm.n_active_params, apm.n_active_params
                ),
                0,
                0,
            )
        else:  # need to set part of the var_cov matrix e.g. if only refined some params
            # first work out the order in self._var_cov
            cumul_pos_dict = {}
            n_cumul_params = 0
            for name, component in six.iteritems(self.components):
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
                    self._var_cov.assign_block(
                        sub, cumul_pos_dict[name], cumul_pos_dict[name2]
                    )

    def update_for_minimisation(self, apm, block_id):
        """Update the scale factors and Ih for the next minimisation iteration."""
        apm_i = apm.apm_list[0]
        basis_fn = self._basis_function.calculate_scales_and_derivatives(
            apm_i, block_id
        )
        self.Ih_table.set_derivatives(basis_fn[1], block_id)
        self.Ih_table.set_inverse_scale_factors(basis_fn[0], block_id)
        self.Ih_table.update_weights(block_id)
        self.Ih_table.calc_Ih(block_id)

    def combine_intensities(self, use_Imid=None):
        """Combine prf and sum intensities to give optimal intensities."""
        try:
            if use_Imid is not None:
                logger.info(
                    "Using previously determined optimal intensity choice: %s\n",
                    OrderedDict(
                        [
                            (use_Imid, str(round(use_Imid, 4))),
                            (0, "profile intensities"),
                            (1, "summation intensities"),
                        ]
                    )[use_Imid],
                )
            else:
                logger.info("Performing profile/summation intensity optimisation.")
            combiner = SingleDatasetIntensityCombiner(self, use_Imid)
        except DialsMergingStatisticsError as e:
            logger.info("Intensity combination failed with the error %s", e)
        else:
            intensity, variance = combiner.calculate_suitable_combined_intensities()
            # update data in reflection table
            self._reflection_table["intensity"].set_selected(
                self.suitable_refl_for_scaling_sel.iselection(), intensity
            )
            self._reflection_table["variance"].set_selected(
                self.suitable_refl_for_scaling_sel.iselection(), variance
            )
            # now set in global_Ih_table
            self.global_Ih_table.update_data_in_blocks(intensity, 0, column="intensity")
            self.global_Ih_table.update_data_in_blocks(variance, 0, column="variance")
            self.global_Ih_table.calc_Ih()
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
                    cpp_calc_sigmasq(jacobian.transpose(), self._var_cov)
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
            logger.info(
                "Scale factors determined during minimisation have now been\n"
                "applied to all reflections for dataset %s.\n",
                self.reflection_table["id"][0],
            )

    def update_error_model(self, error_model, update_Ih=True):
        """Apply a correction to try to improve the error estimate."""
        self.error_model = error_model
        if update_Ih:
            self.global_Ih_table.update_error_model(error_model)
        self.experiment.scaling_model.set_error_model(error_model)

    def adjust_variances(self, caller=None):
        """Apply an aimless-like error model to the variances."""
        error_model = self.experiment.scaling_model.error_model
        if error_model and self.params.weighting.output_optimised_vars:
            # Note : this action has overwritten the variances, so no further
            # error model adjustment should take place, without reinitialising from
            # the input variances (i.e. intensity.prf.variance).
            new_var = error_model.update_variances(
                self.reflection_table["variance"], self.reflection_table["intensity"]
            )
            self.reflection_table["variance"] = new_var
            if caller is None:
                msg = (
                    "The error model has been used to adjust the variances for dataset {0}. \n"
                ).format(self.reflection_table["id"][0])
                logger.info(msg)
        # now increase the errors slightly to take into account the uncertainty in the
        # inverse scale factors
        fractional_error = (
            self.reflection_table["inverse_scale_factor_variance"] ** 0.5
            / self.reflection_table["inverse_scale_factor"]
        )
        variance_scaling = (
            flex.double(self.reflection_table.size(), 1.0) + fractional_error
        )
        self.reflection_table["variance"] *= variance_scaling
        if caller is None:
            msg = (
                "The variances have been adjusted to account for the uncertainty \n"
                "in the scaling model for dataset {0}. \n"
            ).format(self.reflection_table["id"][0])
            logger.info(msg)

    def _select_reflections_for_scaling(self):
        """Select a subset of reflections to use in minimisation."""
        if self.params.reflection_selection.method == "quasi_random":
            block = self.global_Ih_table.Ih_table_blocks[0]
            loc_indices = block.Ih_table["loc_indices"]
            block.Ih_table["s1c"] = (
                self.reflection_table["s1c"]
                .select(self.suitable_refl_for_scaling_sel)
                .select(loc_indices)
            )
            suitable_table = self.reflection_table.select(
                self.suitable_refl_for_scaling_sel
            )
            presel = calculate_scaling_subset_ranges(
                suitable_table, self.params, print_summary=True
            )
            preselection = presel.select(block.Ih_table["loc_indices"])
            self.scaling_selection = calculate_scaling_subset_connected(
                block, self.experiment, self.params, preselection, print_summary=True
            )
        elif self.params.reflection_selection.method == "intensity_ranges":
            overall_scaling_selection = calculate_scaling_subset_ranges_with_E2(
                self.reflection_table, self.params
            )
            self.scaling_selection = overall_scaling_selection.select(
                self.suitable_refl_for_scaling_sel
            )
        elif self.params.reflection_selection.method == "use_all":
            self.scaling_selection = flex.bool(self.n_suitable_refl, True)
        else:
            raise ValueError("Invalid choice for 'reflection_selection.method'.")
        self.scaling_subset_sel = copy.deepcopy(self.scaling_selection)
        self.scaling_selection &= ~self.outliers  # now apply outliers

    def _update_model_data(self):
        """Use the data in the Ih_table to update the model data."""
        assert self.Ih_table is not None
        block_selections = self.Ih_table.get_block_selections_for_dataset(dataset=0)
        for component in self.components.values():
            component.update_reflection_data(block_selections=block_selections)

    def _create_Ih_table(self):
        """Create an Ih_table from the reflection table using the scaling selection."""
        free_set_percentage = 0.0
        if self.params.scaling_options.use_free_set:
            free_set_percentage = self.params.scaling_options.free_set_percentage
        self._Ih_table = IhTable(
            [
                self.reflection_table.select(self.suitable_refl_for_scaling_sel).select(
                    self.scaling_selection
                )
            ],
            self.space_group,
            indices_lists=[self.scaling_selection.iselection()],
            nblocks=self.params.scaling_options.nproc,
            free_set_percentage=free_set_percentage,
            free_set_offset=self.params.scaling_options.free_set_offset,
        )
        if self.error_model:
            variance = self.reflection_table["variance"].select(
                self.suitable_refl_for_scaling_sel
            )
            intensity = self.reflection_table["intensity"].select(
                self.suitable_refl_for_scaling_sel
            )
            new_vars = self.error_model.update_variances(variance, intensity)
            self._Ih_table.update_data_in_blocks(new_vars, 0, column="variance")

    def _configure_model_and_datastructures(self):
        """
        Store the relevant data in the scaling model components.

        This takes the columns from the 'suitable' part of the reflection table (
        which will include outliers). Then a global_Ih_table is created, which can
        be used for outlier rejection. When calculations are done with the model
        components, the correct reflections should first be selected out of the
        stored data.
        """
        sel_reflections = self._reflection_table.select(
            self.suitable_refl_for_scaling_sel
        )
        self.experiment.scaling_model.configure_components(
            sel_reflections, self.experiment, self.params
        )
        self._global_Ih_table = IhTable([sel_reflections], self.space_group, nblocks=1)
        rows = [[key, str(val.n_params)] for key, val in six.iteritems(self.components)]
        st = simple_table(rows, ["correction", "n_parameters"])
        logger.info("The following corrections will be applied to this dataset: \n")
        logger.info(st.format())

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
            self.outliers = flex.bool(self.n_suitable_refl, False)
            self.outliers.set_selected(outlier_indices, True)

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

    def clean_reflection_tables(self):
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

    def set_outliers(self):
        """Set the scaling outliers in the reflection table."""
        outliers = self.outliers
        suitable_isel = self.suitable_refl_for_scaling_sel.iselection()
        outlier_isel = suitable_isel.select(outliers)
        outliers_mask = flex.bool(self.suitable_refl_for_scaling_sel.size(), False)
        outliers_mask.set_selected(outlier_isel, True)
        self.reflection_table.set_flags(
            outliers_mask, self.reflection_table.flags.outlier_in_scaling
        )


class MultiScalerBase(ScalerBase):
    """Base class for scalers handling multiple datasets."""

    def __init__(self, params, experiments, single_scalers):
        """Initialise from a list of single scalers."""
        super(MultiScalerBase, self).__init__()
        self.single_scalers = single_scalers
        self._experiment = experiments[0]
        self._space_group = single_scalers[0].space_group
        self.active_scalers = None
        self._initial_keys = self.single_scalers[0].initial_keys
        self._params = params

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
        if 0 in n_list:
            self._experiment = scalers[0].experiments
        assert len(scalers) == initial_number - len(n_list)
        logger.info("Removed datasets: %s", n_list)

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
        self.global_Ih_table.calc_Ih()
        logger.info(
            "Scale factors determined during minimisation have now been\n"
            "applied to all datasets.\n"
        )

    def adjust_variances(self):
        """Update variances of individual reflection tables."""
        for scaler in self.active_scalers:
            scaler.adjust_variances(caller=self)
        if (
            self.single_scalers[0].experiment.scaling_model.error_model
            and self.params.weighting.output_optimised_vars
        ):
            logger.info(
                "The error model has been used to adjust the variances for all \n"
                "applicable datasets. \n"
            )
        logger.info(
            "The variances have been adjusted to account for the uncertainty \n"
            "in the scaling model for all datasets. \n"
        )

    def clean_reflection_tables(self):
        """Remove unneccesary columns added to reflection tables."""
        for scaler in self.active_scalers:
            scaler.clean_reflection_tables()

    def update_for_minimisation(self, apm, block_id, calc_Ih=True):
        """Update the scale factors and Ih for the next iteration of minimisation."""
        scales = flex.double([])
        derivs = []
        for apm_i in apm.apm_list:
            basis_fn = self._basis_function.calculate_scales_and_derivatives(
                apm_i, block_id
            )
            scales.extend(basis_fn[0])
            derivs.append(basis_fn[1])
        deriv_matrix = sparse.matrix(scales.size(), apm.n_active_params)
        start_row_no = 0
        for j, deriv in enumerate(derivs):
            deriv_matrix.assign_block(deriv, start_row_no, apm.apm_data[j]["start_idx"])
            start_row_no += deriv.n_rows
        self.Ih_table.set_inverse_scale_factors(scales, block_id)
        self.Ih_table.set_derivatives(deriv_matrix, block_id)
        self.Ih_table.update_weights(block_id)
        if calc_Ih:
            self.Ih_table.calc_Ih(block_id)
        # The parallelisation below would work if sparse matrices were
        # pickleable (I think!) - with more benefit for larger number of datasets.'''
        """def task_wrapper(block):
      bf = basis_function(block)
      s, d = bf.calculate_scales_and_derivatives()
      return s, d
    blocks = apm.apm_list
    task_results = easy_mp.parallel_map(func=task_wrapper, iterable=blocks,
      processes=n_datasets, method="multiprocessing",
      preserve_exception_message=True)
    scales_list, derivs_list = zip(*task_results)"""

    def update_error_model(self, error_model, update_Ih=True):
        """Update the error model in Ih table."""
        self.error_model = error_model
        if update_Ih:
            self.global_Ih_table.update_error_model(error_model)
        for scaler in self.active_scalers:
            scaler.experiment.scaling_model.set_error_model(error_model)

    def _update_model_data(self):
        for i, scaler in enumerate(self.active_scalers):
            block_selections = self.Ih_table.get_block_selections_for_dataset(i)
            for component in scaler.components.values():
                component.update_reflection_data(block_selections=block_selections)

    def set_outliers(self):
        """Set the outlier flags in the reflection tables."""
        for scaler in self.active_scalers:
            outliers = scaler.outliers
            suitable_isel = scaler.suitable_refl_for_scaling_sel.iselection()
            outlier_isel = suitable_isel.select(outliers)
            outliers_mask = flex.bool(
                scaler.suitable_refl_for_scaling_sel.size(), False
            )
            outliers_mask.set_selected(outlier_isel, True)
            scaler.reflection_table.set_flags(
                outliers_mask, scaler.reflection_table.flags.outlier_in_scaling
            )

    def _create_global_Ih_table(self):
        tables = [
            s.reflection_table.select(s.suitable_refl_for_scaling_sel)
            for s in self.active_scalers
        ]
        self._global_Ih_table = IhTable(
            tables, self.space_group, nblocks=1, additional_cols=["partiality"]
        )

    def _create_Ih_table(self):
        """Create a new Ih table from the reflection tables."""
        free_set_percentage = 0.0
        if self.params.scaling_options.use_free_set:
            free_set_percentage = self.params.scaling_options.free_set_percentage
        tables = [
            s.reflection_table.select(s.suitable_refl_for_scaling_sel).select(
                s.scaling_selection
            )
            for s in self.active_scalers
        ]
        indices_lists = [s.scaling_selection.iselection() for s in self.active_scalers]
        self._Ih_table = IhTable(
            tables,
            self.space_group,
            indices_lists=indices_lists,
            nblocks=self.params.scaling_options.nproc,
            free_set_percentage=free_set_percentage,
            free_set_offset=self.params.scaling_options.free_set_offset,
        )
        if self.error_model:
            for i, scaler in enumerate(self.active_scalers):
                variance = scaler.reflection_table["variance"].select(
                    scaler.suitable_refl_for_scaling_sel
                )
                intensity = scaler.reflection_table["intensity"].select(
                    scaler.suitable_refl_for_scaling_sel
                )
                new_vars = self.error_model.update_variances(variance, intensity)
                self._Ih_table.update_data_in_blocks(new_vars, i, column="variance")

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
        self._create_Ih_table()
        self._update_model_data()

    @Subject.notify_event(event="performed_outlier_rejection")
    def round_of_outlier_rejection(self, target=None):
        """
        Perform a round of outlier rejection across all datasets.

        After identifying outliers, set the outliers property in individual scalers.
        """
        assert self.active_scalers is not None
        if not self.global_Ih_table:
            self._create_global_Ih_table()
        if self.params.scaling_options.outlier_rejection:
            outlier_index_arrays = determine_outlier_index_arrays(
                self.global_Ih_table,
                self.params.scaling_options.outlier_rejection,
                self.params.scaling_options.outlier_zmax,
                target=target,
            )
            for outlier_indices, scaler in zip(
                outlier_index_arrays, self.active_scalers
            ):
                scaler.outliers = flex.bool(scaler.n_suitable_refl, False)
                scaler.outliers.set_selected(outlier_indices, True)
        logger.debug("Finished outlier rejection.")
        log_memory_usage()

    def _select_reflections_for_scaling(self):
        if self.params.reflection_selection.method == "quasi_random":
            qr = self.params.reflection_selection.quasi_random
            indices, dataset_ids, _ = select_connected_reflections_across_datasets(
                self.global_Ih_table,
                qr.multi_dataset.min_per_dataset,
                qr.multi_dataset.min_multiplicity,
                qr.multi_dataset.Isigma_cutoff,
            )
            header = [
                "Dataset id",
                "reflections \nconnected to \nother datasets",
                "reflections \nhighly connected \nwithin dataset",
                "combined number \nof reflections",
            ]
            rows = []
            for i, scaler in enumerate(self.active_scalers):
                sel = dataset_ids == i
                indices_for_dataset = indices.select(sel)
                scaler.scaling_selection = flex.bool(scaler.n_suitable_refl, False)
                scaler.scaling_selection.set_selected(indices_for_dataset, True)
                # now find good ones from resolution method.
                sel = (
                    self.global_Ih_table.Ih_table_blocks[0].Ih_table["dataset_id"] == i
                )
                indiv_Ih_block = self.global_Ih_table.Ih_table_blocks[0].select(sel)
                loc_indices = indiv_Ih_block.Ih_table["loc_indices"]
                indiv_Ih_block.Ih_table["s1c"] = (
                    scaler.reflection_table["s1c"]
                    .select(scaler.suitable_refl_for_scaling_sel)
                    .select(loc_indices)
                )
                suitable_table = scaler.reflection_table.select(
                    scaler.suitable_refl_for_scaling_sel
                )
                presel = calculate_scaling_subset_ranges(suitable_table, self.params)
                preselection = presel.select(indiv_Ih_block.Ih_table["loc_indices"])

                sel = calculate_scaling_subset_connected(
                    indiv_Ih_block, scaler.experiment, self.params, preselection
                )
                scaler.scaling_selection |= sel
                rows.append(
                    [
                        scaler.experiment.identifier,
                        str(indices_for_dataset.size()),
                        str(sel.count(True)),
                        str(scaler.scaling_selection.count(True)),
                    ]
                )
                scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                scaler.scaling_selection &= ~scaler.outliers
            st = simple_table(rows, header)
            logger.info(
                "Summary of reflections chosen for minimisation from each dataset:"
            )
            logger.info(st.format())
        elif self.params.reflection_selection.method == "intensity_ranges":
            for scaler in self.active_scalers:
                overall_scaling_selection = calculate_scaling_subset_ranges_with_E2(
                    scaler.reflection_table, scaler.params
                )
                scaler.scaling_selection = overall_scaling_selection.select(
                    scaler.suitable_refl_for_scaling_sel
                )
                scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                scaler.scaling_selection &= ~scaler.outliers
        elif self.params.reflection_selection.method == "use_all":
            for scaler in self.active_scalers:
                scaler.scaling_selection = flex.bool(scaler.n_suitable_refl, True)
                scaler.scaling_subset_sel = copy.deepcopy(scaler.scaling_selection)
                scaler.scaling_selection &= ~scaler.outliers
        else:
            raise ValueError("Invalid choice for 'reflection_selection.method'.")


class MultiScaler(MultiScalerBase):
    """Scaler for multiple datasets where all datasets are being minimised."""

    id_ = "multi"

    def __init__(self, params, experiments, single_scalers):
        """
        Initialise a multiscaler from a list of single scalers.

        Create a global_Ih_table, an Ih_table to use for minimisation and update
        the data in the model components.
        """
        logger.info("Configuring a MultiScaler to handle the individual Scalers. \n")
        super(MultiScaler, self).__init__(params, experiments, single_scalers)
        logger.info("Determining symmetry equivalent reflections across datasets.\n")
        self.active_scalers = self.single_scalers
        self._create_global_Ih_table()
        # now select reflections from across the datasets
        self._select_reflections_for_scaling()
        self._create_Ih_table()
        # now add data to scale components from datasets
        self._update_model_data()
        logger.info("Completed configuration of MultiScaler. \n\n" + "=" * 80 + "\n")
        log_memory_usage()

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
                    intensity, variance = combiner.calculate_suitable_combined_intensities(
                        i
                    )
                    scaler.reflection_table["intensity"].set_selected(
                        scaler.suitable_refl_for_scaling_sel.iselection(), intensity
                    )
                    scaler.reflection_table["variance"].set_selected(
                        scaler.suitable_refl_for_scaling_sel.iselection(), variance
                    )
                    self.global_Ih_table.update_data_in_blocks(
                        intensity, i, column="intensity"
                    )
                    self.global_Ih_table.update_data_in_blocks(
                        variance, i, column="variance"
                    )
                    scaler.experiment.scaling_model.record_intensity_combination_Imid(
                        combiner.max_key
                    )
                self.global_Ih_table.calc_Ih()
        else:
            for scaler in self.single_scalers:
                scaler.combine_intensities()


class TargetScaler(MultiScalerBase):
    """A target scaler for scaling datasets against already scaled data."""

    id_ = "target"

    def __init__(self, params, scaled_experiments, scaled_scalers, unscaled_scalers):
        """
        Initialise a multiscaler from a list of single and unscaled scalers.

        First, set the active scalers (the unscaled scalers) and use these to
        create a global_Ih_table. Then, use the scaled_scalers to create a
        target_Ih_table. Create an Ih_table to use for minimisation and use
        the target_Ih_table to set the Ih_values. Finally, update the data in
        the model components.
        """
        logger.info("\nInitialising a TargetScaler instance. \n")
        super(TargetScaler, self).__init__(params, scaled_experiments, scaled_scalers)
        logger.info("Determining symmetry equivalent reflections across datasets.\n")
        self.unscaled_scalers = unscaled_scalers
        self.active_scalers = unscaled_scalers
        self._create_global_Ih_table()
        self._select_reflections_for_scaling()
        tables = [
            s.reflection_table.select(s.suitable_refl_for_scaling_sel).select(
                s.scaling_selection
            )
            for s in self.single_scalers
        ]
        self._target_Ih_table = IhTable(
            tables, self.space_group, nblocks=1
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
        super(TargetScaler, self)._create_Ih_table()
        for block in self._Ih_table.blocked_data_list:
            # this step reduces the number of reflections in each block
            block.match_Ih_values_to_target(self._target_Ih_table)
        self.Ih_table.generate_block_selections()

    def round_of_outlier_rejection(self):
        """Perform a round of targeted outlier rejection."""
        super(TargetScaler, self).round_of_outlier_rejection(
            target=self._target_Ih_table
        )

    def update_for_minimisation(self, apm, block_id, calc_Ih=False):
        """Calcalate the new parameters but don't calculate a new Ih."""
        super(TargetScaler, self).update_for_minimisation(
            apm, block_id, calc_Ih=calc_Ih
        )

    def perform_scaling(
        self, target_type=ScalingTargetFixedIH, engine=None, max_iterations=None
    ):
        """Minimise the scaling model, using a fixed-Ih target."""
        super(TargetScaler, self).perform_scaling(
            target_type=target_type, engine=engine, max_iterations=max_iterations
        )


class NullScaler(ScalerBase):
    """A singlescaler to allow targeted scaling against calculated intensities."""

    id_ = "null"

    def __init__(self, params, experiment, reflection):
        """Set the required properties to use as a scaler for targeted scaling."""
        super(NullScaler, self).__init__()
        self._experiment = experiment
        self._params = params
        self._space_group = self.experiment.crystal.get_space_group()
        self._reflection_table = reflection
        self._initial_keys = list(self._reflection_table.keys())
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
