"""
Definitions of outlier rejection algorithms.

These algorithms use the Ih_table datastructures to perform calculations
in groups of symmetry equivalent reflections. Two functions are provided,
reject_outliers, to reject outlier and set flags given a reflection table
and experiment object, and determine_outlier_index_arrays, which takes an
Ih_table and returns flex.size_t index arrays of the outlier positions.
"""

from __future__ import annotations

import copy
import logging

import numpy as np

from dxtbx import flumpy
from scitbx.array_family import flex

from dials.algorithms.scaling.Ih_table import IhTable
from dials.util.normalisation import quasi_normalisation
from dials_scaling_ext import determine_outlier_indices, limit_outlier_weights

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

    if "inverse_scale_factor" not in reflection_table:
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
        flumpy.from_numpy(outlier_indices), reflection_table.flags.outlier_in_scaling
    )

    return reflection_table


def determine_Esq_outlier_index_arrays(Ih_table, experiment, emax=10.0):
    # first calculate normalised intensities and set in the Ih_table.
    intensities = Ih_table.as_miller_array(experiment.crystal.get_unit_cell())
    normalised_intensities = quasi_normalisation(intensities)

    sel = normalised_intensities.data() > (emax**2)
    n_e2_outliers = sel.count(True)
    if n_e2_outliers:
        logger.info(
            f"{n_e2_outliers} outliers identified from normalised intensity analysis (E\xb2 > {(emax ** 2)})"
        )
    outlier_indices = (
        Ih_table.Ih_table_blocks[0]
        .Ih_table["loc_indices"]
        .iloc[flumpy.to_numpy(sel)]
        .to_numpy()
    )
    if Ih_table.n_datasets == 1:
        return [outlier_indices]
    datasets = (
        Ih_table.Ih_table_blocks[0]
        .Ih_table["dataset_id"]
        .iloc[flumpy.to_numpy(sel)]
        .to_numpy()
    )
    final_outlier_arrays = []
    for i in range(Ih_table.n_datasets):
        final_outlier_arrays.append(outlier_indices[datasets == i])
    return final_outlier_arrays


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
    outlier_rej = None
    if method == "standard":
        outlier_rej = NormDevOutlierRejection(Ih_table, zmax)
    elif method == "simple":
        outlier_rej = SimpleNormDevOutlierRejection(Ih_table, zmax)
    elif method == "target":
        assert target is not None
        outlier_rej = TargetedOutlierRejection(Ih_table, zmax, target)
    elif method is not None:
        raise ValueError(f"Invalid choice of outlier rejection method: {method}")
    if not outlier_rej:
        return [flex.size_t([]) for _ in range(Ih_table.n_datasets)]
    outlier_rej.run()
    outlier_index_arrays = outlier_rej.final_outlier_arrays
    if Ih_table.n_datasets > 1:
        msg = (
            "Combined outlier rejection has been performed across multiple datasets, \n"
        )
    else:
        msg = "A round of outlier rejection has been performed, \n"
    n_outliers = sum(len(i) for i in outlier_index_arrays)
    msg += f"{n_outliers} outliers have been identified. \n"
    logger.info(msg)
    return outlier_index_arrays


class OutlierRejectionBase:
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
        self._datasets = np.array([], dtype=np.uint64).reshape((0,))  # flex.int([])
        self._zmax = zmax
        self._outlier_indices = np.array([], dtype=np.uint64).reshape(
            (0,)
        )  # flex.size_t([])
        self.final_outlier_arrays = None

    def run(self):
        """Run the outlier rejection algorithm, implemented by a subclass."""
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
            return [self._outlier_indices]
        final_outlier_arrays = []
        for i in range(self._n_datasets):
            final_outlier_arrays.append(self._outlier_indices[self._datasets == i])
        return final_outlier_arrays

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        raise NotImplementedError()


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
        super().__init__(Ih_table, zmax)

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        Ih_table = self._Ih_table_block
        target = self._target_Ih_table_block
        target_asu_Ih_dict = dict(
            zip(target.asu_miller_index, zip(target.Ih_values, target.variances))
        )
        target_Ih_value = np.zeros(Ih_table.size)
        target_Ih_sigmasq = np.zeros(Ih_table.size)
        for j, miller_idx in enumerate(Ih_table.asu_miller_index):
            try:
                vals = target_asu_Ih_dict[miller_idx]
            except KeyError:
                pass
            else:
                target_Ih_value[j] = vals[0]
                target_Ih_sigmasq[j] = vals[1]

        nz_sel = target_Ih_value != 0.0
        target_Ih_value = target_Ih_value[nz_sel]
        target_Ih_sigmasq = target_Ih_sigmasq[nz_sel]
        Ih_table = Ih_table.select(nz_sel)
        norm_dev = (
            Ih_table.intensities - (Ih_table.inverse_scale_factors * target_Ih_value)
        ) / (
            np.sqrt(
                Ih_table.variances
                + (np.square(Ih_table.inverse_scale_factors) * target_Ih_sigmasq)
            )
        )
        outliers_sel = np.abs(norm_dev) > self._zmax
        outliers_isel = np.nonzero(nz_sel)[0][outliers_sel]

        outliers = np.full(self._Ih_table_block.size, False)
        outliers[outliers_isel] = True

        self._outlier_indices = np.concatenate(
            [
                self._outlier_indices,
                self._Ih_table_block.Ih_table["loc_indices"].iloc[outliers],
            ]
        )
        self._datasets = np.concatenate(
            [
                self._datasets,
                self._Ih_table_block.Ih_table["dataset_id"].iloc[outliers],
            ]
        )


class SimpleNormDevOutlierRejection(OutlierRejectionBase):
    """Algorithm using normalised deviations from the weighted intensity means.

    In this case, the weighted mean is calculated from all reflections in
    the symmetry group excluding the test reflection.
    """

    def __init__(self, Ih_table, zmax):
        super().__init__(Ih_table, zmax)
        self.weights = flumpy.to_numpy(
            limit_outlier_weights(
                copy.deepcopy(self._Ih_table_block.weights),
                self._Ih_table_block.h_index_matrix,
            )
        )

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        Ih_table = self._Ih_table_block
        intensity = Ih_table.intensities
        g = Ih_table.inverse_scale_factors
        w = self.weights
        wgIsum = Ih_table.sum_in_groups(w * g * intensity, output="per_refl")
        wg2sum = Ih_table.sum_in_groups(w * g * g, output="per_refl")

        # guard against zero division errors - can happen due to rounding errors
        # or bad data giving g values are very small
        zero_sel = wg2sum == 0.0
        # set as one for now, then mark as outlier below. This will only affect if
        # g is near zero, if w is zero then throw an assertionerror.
        wg2sum[zero_sel] = 1.0

        assert np.all(w > 0)  # guard against division by zero
        norm_dev = (intensity - (g * wgIsum / wg2sum)) / (
            np.sqrt((1.0 / w) + np.square(g / wg2sum))
        )
        norm_dev[zero_sel] = 1000  # to trigger rejection
        outliers = np.abs(norm_dev) > self._zmax

        self._outlier_indices = np.concatenate(
            [
                self._outlier_indices,
                self._Ih_table_block.Ih_table["loc_indices"].iloc[outliers].to_numpy(),
            ]
        )
        self._datasets = np.concatenate(
            [
                self._datasets,
                self._Ih_table_block.Ih_table["dataset_id"].iloc[outliers].to_numpy(),
            ]
        )


class NormDevOutlierRejection(OutlierRejectionBase):
    """Algorithm using normalised deviations from the weighted intensity means.

    In this case, the weighted mean is calculated from all reflections in
    the symmetry group excluding the test reflection.
    """

    def __init__(self, Ih_table, zmax):
        super().__init__(Ih_table, zmax)
        self.weights = flumpy.to_numpy(
            limit_outlier_weights(
                copy.deepcopy(self._Ih_table_block.weights),
                self._Ih_table_block.h_index_matrix,
            )
        )

    def _do_outlier_rejection(self):
        """Add indices (w.r.t. the Ih_table data) to self._outlier_indices."""
        self._round_of_outlier_rejection()
        n_outliers = self._outlier_indices.size
        n_new_outliers = n_outliers
        while n_new_outliers:
            self._round_of_outlier_rejection()
            n_new_outliers = self._outlier_indices.size - n_outliers
            n_outliers = self._outlier_indices.size

    def _round_of_outlier_rejection(self):
        """
        Calculate normal deviations from the data in the Ih_table.
        """
        Ih_table = self._Ih_table_block
        intensity = Ih_table.intensities
        g = Ih_table.inverse_scale_factors
        w = self.weights
        wgIsum = Ih_table.sum_in_groups(w * g * intensity, output="per_refl")
        wg2sum = Ih_table.sum_in_groups(w * g * g, output="per_refl")
        wgIsum_others = wgIsum - (w * g * intensity)
        wg2sum_others = wg2sum - (w * g * g)
        # Now do the rejection analysis if n_in_group > 2
        nh = Ih_table.calc_nh()
        sel = nh > 2
        wg2sum_others_sel = wg2sum_others[sel]
        wgIsum_others_sel = wgIsum_others[sel]

        # guard against zero division errors - can happen due to rounding errors
        # or bad data giving g values are very small
        zero_sel = wg2sum_others_sel == 0.0
        # set as one for now, then mark as outlier below. This will only affect if
        # g is near zero, if w is zero then throw an assertionerror.
        wg2sum_others_sel[zero_sel] = 1.0
        g_sel = g[sel]
        I_sel = intensity[sel]
        w_sel = w[sel]

        assert np.all(w_sel > 0)  # guard against division by zero
        norm_dev = (I_sel - (g_sel * wgIsum_others_sel / wg2sum_others_sel)) / (
            np.sqrt((1.0 / w_sel) + (np.square(g_sel) / wg2sum_others_sel))
        )
        norm_dev[zero_sel] = 1000  # to trigger rejection
        z_score = np.abs(norm_dev)
        # Want an array same size as Ih table.
        all_z_scores = np.zeros(Ih_table.size)
        # all_z_scores.set_selected(sel.iselection(), z_score)
        all_z_scores[sel] = z_score
        outlier_indices, other_potential_outliers = determine_outlier_indices(
            Ih_table.h_index_matrix, flumpy.from_numpy(all_z_scores), self._zmax
        )
        sel = np.full(Ih_table.size, False, dtype=bool)
        outlier_indices = flumpy.to_numpy(outlier_indices)
        sel[outlier_indices] = True
        lsel = self._Ih_table_block.Ih_table["loc_indices"].iloc[sel].to_numpy()
        dsel = self._Ih_table_block.Ih_table["dataset_id"].iloc[sel].to_numpy()
        self._outlier_indices = np.concatenate(
            [
                self._outlier_indices,
                lsel,
            ]
        )
        self._datasets = np.concatenate(
            [
                self._datasets,
                dsel,
            ]
        )

        sel = np.full(Ih_table.size, False, dtype=bool)
        sel[flumpy.to_numpy(other_potential_outliers)] = True
        self._Ih_table_block = self._Ih_table_block.select(sel)
        self.weights = self.weights[sel]
