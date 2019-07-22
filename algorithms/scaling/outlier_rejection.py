"""
Definitions of outlier rejection algorithms.

These algorithms use the Ih_table datastructures to perform calculations
in groups of symmetry equivalent reflections. Two functions are provided,
reject_outliers, to reject outlier and set flags given a reflection table
and experiment object, and determine_outlier_index_arrays, which takes an
Ih_table and returns flex.size_t index arrays of the outlier positions.
"""
from __future__ import absolute_import, division, print_function
import abc
import logging
from scitbx.array_family import flex
from dials.algorithms.scaling.Ih_table import IhTable
from dials_scaling_ext import determine_outlier_indices

logger = logging.getLogger("dials")


def reject_outliers(reflection_table, experiment, method="standard", zmax=6.0):
    """
    Run an outlier algorithm on symmetry-equivalent intensities.

    This method runs an intensity-based outlier rejection algorithm, comparing
    the deviations from the weighted mean in groups of symmetry equivalent
    reflections. The outliers are determined and the outlier_in_scaling flag
    is set in the reflection table.

    The values intensity and variance must be set in the reflection table;
    these should be corrected but unscaled values, as an inverse_scale_factor
    will be applied during outlier rejection if this is present in the reflection
    table. The reflection table should also be prefiltered (e.g. not-integrated
    reflections should not be present) as no further filtering is done on the
    input table.

    Args:
        reflection_table: A reflection table.
        experiment: A single experiment object.
        method (str): Name (alias) of outlier rejection algorithm to use.
        zmax (float): Normalised deviation threshold for classifying an outlier.

    Returns:
        reflection_table: The input table with the outlier_in_scaling flag set.

    """
    assert "intensity" in reflection_table, "reflection table has no 'intensity' column"
    assert "variance" in reflection_table, "reflection table has no 'variance' column"

    if not "inverse_scale_factor" in reflection_table:
        reflection_table["inverse_scale_factor"] = flex.double(
            reflection_table.size(), 1.0
        )

    Ih_table = IhTable(
        [reflection_table], experiment.crystal.get_space_group(), nblocks=1
    )
    outlier_indices = determine_outlier_index_arrays(
        Ih_table, method=method, zmax=zmax
    )[0]

    # Unset any existing outlier flags before setting the new ones
    reflection_table.unset_flags(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling),
        reflection_table.flags.outlier_in_scaling,
    )
    reflection_table.set_flags(
        outlier_indices, reflection_table.flags.outlier_in_scaling
    )

    return reflection_table


def determine_outlier_index_arrays(Ih_table, method="standard", zmax=6.0, target=None):
    """
    Run an outlier algorithm and return the outlier indices.

    Args:
        Ih_table: A dials.algorithms.scaling.Ih_table.IhTable.
        method (str): Name (alias) of outlier rejection algorithm to use. If
            method=target, then the optional argument target must also
            be specified. Implemented methods; standard, simple, target.
        zmax (float): Normalised deviation threshold for classifying an outlier.
        target (Optional[IhTable]): An IhTable to use to obtain target Ih for
            outlier rejectiob, if method=target.

    Returns:
        outlier_index_arrays (list): A list of flex.size_t arrays, with one
            array per dataset that was used to create the Ih_table. Importantly,
            the indices are the indices of the reflections in the initial
            reflection table used to create the Ih_table, not the indices of the
            data in the Ih_table.

    Raises:
        ValueError: if an invalid choice is made for the method.

    """
    if method == "standard":
        outlier_index_arrays = NormDevOutlierRejection(
            Ih_table, zmax
        ).final_outlier_arrays
    elif method == "simple":
        outlier_index_arrays = SimpleNormDevOutlierRejection(
            Ih_table, zmax
        ).final_outlier_arrays
    elif method == "target":
        assert target is not None
        outlier_index_arrays = TargetedOutlierRejection(
            Ih_table, zmax, target
        ).final_outlier_arrays
    elif method is None:
        return [flex.size_t([]) for _ in range(Ih_table.n_datasets)]
    else:
        raise ValueError("Invalid choice of outlier rejection method: %s" % method)
    if Ih_table.n_datasets > 1:
        msg = (
            "Combined outlier rejection has been performed across multiple datasets, \n"
        )
    else:
        msg = "A round of outlier rejection has been performed, \n"
    n_outliers = sum([len(i) for i in outlier_index_arrays])
    msg += "{} outliers have been identified. \n".format(n_outliers)
    logger.info(msg)
    return outlier_index_arrays


class OutlierRejectionBase(object):
    """
    Base class for outlier rejection algorithms using an IhTable datastructure.

    Subclasses must implement the _do_outlier_rejection method, which must
    add the indices of outliers to the _outlier_indices attribute. The algorithms
    are run upon initialisation and result in the population of the
    :obj:`final_outlier_arrays`.

    Attributes:
        final_outlier_arrays (:obj:`list`): A list of flex.size_t arrays of outlier
            indices w.r.t. the order of the initial reflection tables used to
            create the Ih_table.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, Ih_table, zmax):
        """Set up and run the outlier rejection algorithm."""
        assert (
            Ih_table.n_work_blocks == 1
        ), """
Outlier rejection algorithms require an Ih_table with nblocks = 1"""
        # Note: could be possible to code for nblocks > 1
        self._Ih_table_block = Ih_table.blocked_data_list[0]
        self._n_datasets = Ih_table.n_datasets
        self._block_selections = Ih_table.blocked_selection_list[0]
        self._ids = self._Ih_table_block.Ih_table["dataset_id"]
        self._zmax = zmax
        self._outlier_indices = flex.size_t([])
        self._do_outlier_rejection()
        self.final_outlier_arrays = self._determine_outlier_indices()

    def _determine_outlier_indices(self):
        """
        Determine outlier indices with respect to the input reflection tables.

        Transform the outlier indices w.r.t the Ih_table, determined during the
        algorithm, to outlier indices w.r.t the initial reflection tables used
        to create the Ih_table, separated by reflection table.

        Returns:
            final_outlier_arrays (:obj:`list`): A list of flex.size_t arrays of
                outlier indices w.r.t. the order of the data in the initial
                reflection tables used to create the Ih_table.

        """
        if self._n_datasets == 1:
            return [self._block_selections[0].select(self._outlier_indices)]
        final_outlier_arrays = []
        ids = self._ids.select(self._outlier_indices)
        offset = 0
        for i in range(self._n_datasets):
            outlier_array_i = self._outlier_indices.select(ids == i) - offset
            final_outlier_arrays.append(
                self._block_selections[i].select(outlier_array_i)
            )
            offset += self._block_selections[i].size()
        return final_outlier_arrays

    @abc.abstractmethod
    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""


class TargetedOutlierRejection(OutlierRejectionBase):
    """Implementation of an outlier rejection algorithm against a target.

    This algorithm requires a target Ih_table in addition to an Ih_table
    for the dataset under investigation. Normalised deviations are
    calculated from the intensity values in the target table.
    """

    def __init__(self, Ih_table, zmax, target):
        """Set a target Ih_table and run the outlier rejection."""
        assert (
            target.n_work_blocks == 1
        ), """
Targeted outlier rejection requires a target Ih_table with nblocks = 1"""
        self._target_Ih_table_block = target.blocked_data_list[0]
        self._target_Ih_table_block.calc_Ih()
        super(TargetedOutlierRejection, self).__init__(Ih_table, zmax)

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        Ih_table = self._Ih_table_block
        target = self._target_Ih_table_block
        target_asu_Ih_dict = dict(
            zip(target.asu_miller_index, zip(target.Ih_values, target.variances))
        )
        Ih_table.Ih_table["target_Ih_value"] = flex.double(Ih_table.size, 0.0)
        Ih_table.Ih_table["target_Ih_sigmasq"] = flex.double(Ih_table.size, 0.0)
        for j, miller_idx in enumerate(Ih_table.asu_miller_index):
            if miller_idx in target_asu_Ih_dict:
                Ih_table.Ih_table["target_Ih_value"][j] = target_asu_Ih_dict[
                    miller_idx
                ][0]
                Ih_table.Ih_table["target_Ih_sigmasq"][j] = target_asu_Ih_dict[
                    miller_idx
                ][1]

        nz_sel = Ih_table.Ih_table["target_Ih_value"] != 0.0
        Ih_table = Ih_table.select(nz_sel)
        norm_dev = (
            Ih_table.intensities
            - (Ih_table.inverse_scale_factors * Ih_table.Ih_table["target_Ih_value"])
        ) / (
            (
                Ih_table.variances
                + (
                    (Ih_table.inverse_scale_factors ** 2)
                    * Ih_table.Ih_table["target_Ih_sigmasq"]
                )
            )
            ** 0.5
        )
        outliers_sel = flex.abs(norm_dev) > self._zmax
        outliers_isel = nz_sel.iselection().select(outliers_sel)
        self._outlier_indices.extend(outliers_isel)


class SimpleNormDevOutlierRejection(OutlierRejectionBase):
    """Algorithm using normalised deviations from the weighted intensity means.

    In this case, the weighted mean is calculated from all reflections in
    the symmetry group excluding the test reflection.
    """

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        Ih_table = self._Ih_table_block
        I = Ih_table.intensities
        g = Ih_table.inverse_scale_factors
        w = Ih_table.weights
        wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
        wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix

        # guard against zero divison errors - can happen due to rounding errors
        # or bad data giving g values are very small
        zero_sel = wg2sum == 0.0
        # set as one for now, then mark as outlier below. This will only affect if
        # g is near zero, if w is zero then throw an assertionerror.
        wg2sum.set_selected(zero_sel, 1.0)

        assert w.all_gt(0)  # guard against division by zero
        norm_dev = (I - (g * wgIsum / wg2sum)) / (
            ((1.0 / w) + ((g / wg2sum) ** 2)) ** 0.5
        )
        norm_dev.set_selected(zero_sel, 1000)  # to trigger rejection
        outliers_sel = flex.abs(norm_dev) > self._zmax

        self._outlier_indices.extend(outliers_sel.iselection())


class NormDevOutlierRejection(OutlierRejectionBase):
    """Algorithm using normalised deviations from the weighted intensity means.

    In this case, the weighted mean is calculated from all reflections in
    the symmetry group excluding the test reflection.
    """

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        outlier_indices, other_potential_outliers = self._round_of_outlier_rejection()
        self._outlier_indices.extend(outlier_indices)
        internal_potential_outliers = other_potential_outliers
        while other_potential_outliers:
            good_sel = flex.bool(self._Ih_table_block.Ih_table.size(), False)
            good_sel.set_selected(internal_potential_outliers, True)
            self._Ih_table_block = self._Ih_table_block.select(good_sel)
            other_potential_outliers, internal_potential_outliers = self._check_for_more_outliers(
                other_potential_outliers
            )

    def _check_for_more_outliers(self, other_potential_outliers):
        """
        Recursive check for further outliers.

        Each iteration creates a new reduced-size Ih_table_block, which retains
        only symmetry groups that need further testing. Outlier indices must be
        transformed to give indices with respect to the initial Ih_table_block.

        Args:
            other_potential_outliers: A flex.size_t array of indices with respect
                to the initial Ih_table data

        """
        # Find outlier indices with respect to reduced Ih_table block
        internal_outlier_indices, internal_other_potential_outliers = (
            self._round_of_outlier_rejection()
        )
        outliers_wrt_original = other_potential_outliers.select(
            internal_outlier_indices
        )
        self._outlier_indices.extend(outliers_wrt_original)
        new_other_potential_outliers = other_potential_outliers.select(
            internal_other_potential_outliers
        )  # still wrt original Ih_table data
        return new_other_potential_outliers, internal_other_potential_outliers

    def _round_of_outlier_rejection(self):
        """
        Calculate normal deviations from the data in the Ih_table.

        Returns:
            (tuple): tuple containing:
                outlier_indices: A flex.size_t array of outlier indices w.r.t
                    the current Ih_table
                other_potential_outliers: A flex.size_t array of indices from
                    the symmetry groups where outliers were found, excluding the
                    indices of the outliers themselves (indices w.r.t current
                    Ih_table).

        """
        Ih_table = self._Ih_table_block
        I = Ih_table.intensities
        g = Ih_table.inverse_scale_factors
        w = Ih_table.weights
        wgIsum = ((w * g * I) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
        wg2sum = ((w * g * g) * Ih_table.h_index_matrix) * Ih_table.h_expand_matrix
        wgIsum_others = wgIsum - (w * g * I)
        wg2sum_others = wg2sum - (w * g * g)
        # Now do the rejection analyis if n_in_group > 2
        nh = Ih_table.calc_nh()
        sel = nh > 2
        wg2sum_others_sel = wg2sum_others.select(sel)
        wgIsum_others_sel = wgIsum_others.select(sel)

        # guard against zero divison errors - can happen due to rounding errors
        # or bad data giving g values are very small
        zero_sel = wg2sum_others_sel == 0.0
        # set as one for now, then mark as outlier below. This will only affect if
        # g is near zero, if w is zero then throw an assertionerror.
        wg2sum_others_sel.set_selected(zero_sel, 1.0)
        g_sel = g.select(sel)
        I_sel = I.select(sel)
        w_sel = w.select(sel)

        assert w_sel.all_gt(0)  # guard against division by zero
        norm_dev = (I_sel - (g_sel * wgIsum_others_sel / wg2sum_others_sel)) / (
            ((1.0 / w_sel) + (g_sel ** 2 / wg2sum_others_sel)) ** 0.5
        )
        norm_dev.set_selected(zero_sel, 1000)  # to trigger rejection
        z_score = flex.abs(norm_dev)
        # Want an array same size as Ih table.
        all_z_scores = flex.double(Ih_table.size, 0.0)
        all_z_scores.set_selected(sel.iselection(), z_score)
        outlier_indices, other_potential_outliers = determine_outlier_indices(
            Ih_table.h_index_matrix, all_z_scores, self._zmax
        )
        return outlier_indices, other_potential_outliers
