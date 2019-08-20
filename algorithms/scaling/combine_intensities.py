"""
Optimise the combination of profile and summation intensity values.
"""
from __future__ import absolute_import, division, print_function
import logging
import boost.python
from libtbx.table_utils import simple_table
from cctbx import miller, crystal
from dials.array_family import flex
from dials.algorithms.scaling.scaling_utilities import DialsMergingStatisticsError

miller_ext = boost.python.import_ext("cctbx_miller_ext")
logger = logging.getLogger("dials")


def fast_merging_stats(array):
    """
    Quickly calculate required merging stats for intensity combination.

    This is a cut-down version of iobtx.merging_statistics.merging_stats.
    """
    assert array.sigmas() is not None
    positive_sel = array.sigmas() > 0
    array = array.select(positive_sel)
    array = array.sort("packed_indices")
    merge_ext = miller_ext.merge_equivalents_obs(
        array.indices(), array.data(), array.sigmas(), use_internal_variance=True
    )
    r_meas = merge_ext.r_meas
    cc_one_half = miller.compute_cc_one_half(unmerged=array, return_n_refl=False)
    return r_meas, cc_one_half


def map_indices_to_asu(miller_indices, space_group):
    """Map the indices to the asymmetric unit."""
    crystal_symmetry = crystal.symmetry(space_group=space_group)
    miller_set = miller.set(
        crystal_symmetry=crystal_symmetry, indices=miller_indices, anomalous_flag=False
    )
    miller_set_in_asu = miller_set.map_to_asu()
    return miller_set_in_asu.indices()


def _make_reflection_table_from_scaler(scaler):
    """Copy across required columns and filter data."""
    reflections = flex.reflection_table()
    required_cols = [
        "intensity.prf.value",
        "intensity.prf.variance",
        "intensity.sum.value",
        "intensity.sum.variance",
        "prescaling_correction",
        "inverse_scale_factor",
        "miller_index",
    ]
    optional_cols = ["partiality"]
    for col in required_cols:
        reflections[col] = scaler.reflection_table[col]
    for col in optional_cols:
        if col in scaler.reflection_table:
            reflections[col] = scaler.reflection_table[col]
    # now select good data
    sel = _get_filter_selection(scaler.reflection_table)
    suitable_isel = scaler.suitable_refl_for_scaling_sel.iselection()
    outlier_isel = suitable_isel.select(scaler.outliers)
    not_outliers = flex.bool(reflections.size(), True)
    not_outliers.set_selected(outlier_isel, False)
    reflections = reflections.select(sel & not_outliers)
    reflections["miller_index"] = map_indices_to_asu(
        reflections["miller_index"], scaler.space_group
    )
    logger.debug("Reflection table size for combining: %s", reflections.size())
    return reflections


class SingleDatasetIntensityCombiner(object):

    """
    Class to combine profile and summation intensities for a single datset.
    """

    def __init__(self, scaler, use_Imid=None):
        if "intensity.prf.value" not in scaler.reflection_table:
            self.max_key = 1
            logger.info(
                "No profile intensities found, skipping profile/summation intensity combination."
            )
            return
        self.scaler = scaler
        self.experiment = scaler.experiment
        if use_Imid is not None:
            self.max_key = use_Imid
        else:
            self.Imids = scaler.params.reflection_selection.combine.Imid
            self.dataset = _make_reflection_table_from_scaler(self.scaler)
            if "partiality" in self.dataset:
                raw_intensities = (
                    self.dataset["intensity.sum.value"].as_double()
                    / self.dataset["partiality"]
                )
            else:
                raw_intensities = self.dataset["intensity.sum.value"].as_double()
            logger.debug("length of raw intensity array: %s", raw_intensities.size())
            self._determine_Imids(raw_intensities)
            header = ["Combination", "CC1/2", "Rmeas"]
            rows, results = self._test_Imid_combinations()
            st = simple_table(rows, header)
            logger.info(st.format())

            self.max_key = min(results, key=results.get)
            if self.max_key == 0:
                logger.info("Profile intensities determined to be best for scaling. \n")
            elif self.max_key == 1:
                logger.info(
                    "Summation intensities determined to be best for scaling. \n"
                )
            else:
                logger.info(
                    "Combined intensities with Imid = %s determined to be best for scaling. \n",
                    self.max_key,
                )

    def calculate_suitable_combined_intensities(self):
        """Combine the 'suitable for scaling' intensities in the scaler."""
        return _calculate_suitable_combined_intensities(self.scaler, self.max_key)

    def _determine_Imids(self, raw_intensities):
        if not self.Imids:
            avg = flex.mean(raw_intensities)
            Imid = flex.max(raw_intensities) / 10.0
            Imid_list = [0, 1, avg, Imid]
            while Imid > avg:
                Imid /= 10.0
                Imid_list.append(Imid)
            self.Imids = Imid_list

    def _test_Imid_combinations(self):
        """Test the different combinations, returning the rows and results dict."""
        rows = []
        results = {}

        for Imid in self.Imids:
            Int, Var = _get_Is_from_Imidval(self.dataset, Imid)
            miller_set = miller.set(
                crystal_symmetry=self.experiment.crystal.get_crystal_symmetry(
                    assert_is_compatible_unit_cell=False
                ),
                indices=self.dataset["miller_index"],
                anomalous_flag=False,
            )
            i_obs = miller.array(
                miller_set,
                data=(
                    Int
                    * self.dataset["prescaling_correction"]
                    / self.dataset["inverse_scale_factor"]
                ),
            )
            i_obs.set_observation_type_xray_intensity()
            i_obs.set_sigmas(
                (Var ** 0.5)
                * self.dataset["prescaling_correction"]
                / self.dataset["inverse_scale_factor"]
            )
            try:
                rmeas, cchalf = fast_merging_stats(array=i_obs)
                logger.debug("Imid: %s, Rmeas %s, cchalf %s", Imid, rmeas, cchalf)
            except RuntimeError:
                raise DialsMergingStatisticsError(
                    "Unable to merge for intensity combination"
                )

            # record the results
            results[Imid] = rmeas
            res_str = {0: "prf only", 1: "sum only"}
            if not Imid in res_str:
                res_str[Imid] = "Imid = " + str(round(Imid, 2))
            rows.append([res_str[Imid], str(round(cchalf, 5)), str(round(rmeas, 5))])

        return rows, results


def _calculate_suitable_combined_intensities(scaler, max_key):
    reflections = scaler.reflection_table
    suitable = scaler.suitable_refl_for_scaling_sel
    suitable_conv = reflections["prescaling_correction"].select(suitable)
    Ipr = reflections["intensity.prf.value"].select(suitable)
    Vpr = reflections["intensity.prf.variance"].select(suitable)
    Isum = reflections["intensity.sum.value"].select(suitable)
    Vsum = reflections["intensity.sum.variance"].select(suitable)
    if "partiality" in reflections:
        inv_p = _determine_inverse_partiality(reflections)
        inv_p = inv_p.select(suitable)
    if max_key == 0:
        intensity = Ipr * suitable_conv
        variance = Vpr * suitable_conv * suitable_conv
    elif max_key == 1:
        if "partiality" in reflections:
            intensity = Isum * suitable_conv * inv_p
            variance = Vsum * (suitable_conv * inv_p) ** 2
        else:
            intensity = Isum * suitable_conv
            variance = Vsum * suitable_conv * suitable_conv
    else:
        if "partiality" in reflections:
            Int, Var = _calculate_combined_raw_intensities(
                Ipr, Isum * inv_p, Vpr, Vsum * inv_p * inv_p, max_key
            )
        else:
            Int, Var = _calculate_combined_raw_intensities(
                Ipr, Isum, Vpr, Vsum, max_key
            )
        intensity = Int * suitable_conv
        variance = Var * suitable_conv * suitable_conv
    return intensity, variance


class MultiDatasetIntensityCombiner(object):

    """
    Class to combine profile and summation intensities for multiple datasets.
    """

    def __init__(self, multiscaler):
        self.active_scalers = multiscaler.active_scalers
        self.experiment = multiscaler.experiment
        self.Imids = multiscaler.params.reflection_selection.combine.Imid
        # first copy across relevant data that's needed
        self.good_datasets = []
        for i, scaler in enumerate(self.active_scalers):
            if "intensity.prf.value" in scaler.reflection_table:
                self.good_datasets.append(i)
        self.datasets = [
            _make_reflection_table_from_scaler(self.active_scalers[i])
            for i in self.good_datasets
        ]
        raw_intensities = self._get_raw_intensity_array()
        logger.debug("length of raw intensity array: %s", raw_intensities.size())
        self._determine_Imids(raw_intensities)

        header = ["Combination", "CC1/2", "Rmeas"]
        rows, results = self._test_Imid_combinations()
        st = simple_table(rows, header)
        logger.info(st.format())

        self.max_key = min(results, key=results.get)
        while results[self.max_key] < 0:
            del results[self.max_key]
            self.max_key = min(results, key=results.get)
        if self.max_key == 0:
            logger.info("Profile intensities determined to be best for scaling. \n")
        elif self.max_key == 1:
            logger.info("Summation intensities determined to be best for scaling. \n")
        else:
            logger.info(
                "Combined intensities with Imid = %s determined to be best for scaling. \n",
                self.max_key,
            )

    def calculate_suitable_combined_intensities(self, dataset):
        """Combine the 'suitable for scaling' intensities in the scaler."""
        if dataset not in self.good_datasets:
            return _calculate_suitable_combined_intensities(
                self.active_scalers[dataset], 1
            )
        return _calculate_suitable_combined_intensities(
            self.active_scalers[dataset], self.max_key
        )

    def _get_raw_intensity_array(self):
        intensities = flex.double()
        for dataset in self.datasets:
            if "partiality" in dataset:
                intensities.extend(
                    dataset["intensity.sum.value"].as_double() / dataset["partiality"]
                )
            else:
                intensities.extend(dataset["intensity.sum.value"].as_double())
        return intensities

    def _determine_Imids(self, raw_intensities):
        if not self.Imids:
            avg = flex.mean(raw_intensities)
            Imid = flex.max(raw_intensities) / 10.0
            Imid_list = [0, 1, avg, Imid]
            while Imid > avg:
                Imid /= 10.0
                Imid_list.append(Imid)
            self.Imids = Imid_list

    def _test_Imid_combinations(self):
        rows = []
        results = {}
        for Imid in self.Imids:
            combined_intensities = flex.double([])
            combined_sigmas = flex.double([])
            combined_scales = flex.double([])
            combined_indices = flex.miller_index([])
            for dataset in self.datasets:
                Int, Var = _get_Is_from_Imidval(dataset, Imid)
                Int *= dataset["prescaling_correction"]
                sigma = (Var ** 0.5) * dataset["prescaling_correction"]
                combined_intensities.extend(Int)
                combined_sigmas.extend(sigma)
                combined_scales.extend(dataset["inverse_scale_factor"])
                combined_indices.extend(dataset["miller_index"])
            # apply scale factor before determining merging stats
            miller_set = miller.set(
                crystal_symmetry=self.experiment.crystal.get_crystal_symmetry(),
                indices=combined_indices,
                anomalous_flag=False,
            )
            i_obs = miller.array(
                miller_set, data=combined_intensities / combined_scales
            )
            i_obs.set_observation_type_xray_intensity()
            i_obs.set_sigmas(combined_sigmas / combined_scales)
            try:
                rmeas, cchalf = fast_merging_stats(array=i_obs)
                logger.debug("Imid: %s, Rmeas %s, cchalf %s", Imid, rmeas, cchalf)
            except RuntimeError:
                raise DialsMergingStatisticsError(
                    "Unable to merge for intensity combination"
                )

            # record the results
            results[Imid] = rmeas
            res_str = {0: "prf only", 1: "sum only"}
            if not Imid in res_str:
                res_str[Imid] = "Imid = " + str(round(Imid, 2))
            rows.append([res_str[Imid], str(round(cchalf, 5)), str(round(rmeas, 5))])
        return rows, results


### Helper functions for combine_intensities


def _get_Is_from_Imidval(reflections, Imid):
    """Intepret the Imid value to extract and return the Icomb and Vcomb values."""
    if Imid == 0:  # special value to trigger prf
        Int = reflections["intensity.prf.value"]
        Var = reflections["intensity.prf.variance"]
    elif Imid == 1:  # special value to trigger sum
        if "partiality" in reflections:
            Int = reflections["intensity.sum.value"] / reflections["partiality"]
            Var = reflections["intensity.sum.variance"] / (
                reflections["partiality"] ** 2
            )
        else:
            Int = reflections["intensity.sum.value"]
            Var = reflections["intensity.sum.variance"]
    else:
        if "partiality" in reflections:
            Int, Var = _calculate_combined_raw_intensities(
                reflections["intensity.prf.value"],
                reflections["intensity.sum.value"] / reflections["partiality"],
                reflections["intensity.prf.variance"],
                reflections["intensity.sum.variance"]
                / (reflections["partiality"] ** 2),
                Imid,
            )
        else:
            Int, Var = _calculate_combined_raw_intensities(
                reflections["intensity.prf.value"],
                reflections["intensity.sum.value"],
                reflections["intensity.prf.variance"],
                reflections["intensity.sum.variance"],
                Imid,
            )
    return Int, Var


def _get_filter_selection(reflections):
    bad_sel = (
        reflections.get_flags(reflections.flags.bad_for_scaling, all=False)
        | (reflections["intensity.prf.variance"] <= 0)
        | (reflections["intensity.sum.variance"] <= 0)
        | (reflections["inverse_scale_factor"] <= 0)
    )
    integrated = reflections.get_flags(reflections.flags.integrated, all=True)
    good_sel = integrated & ~bad_sel
    if "partiality" in reflections:
        good_sel &= reflections["partiality"] > 0
    return good_sel


def _determine_inverse_partiality(reflections):
    inverse_partiality = flex.double(reflections.size(), 1.0)
    nonzero_partiality_sel = reflections["partiality"] > 0.0
    good_refl = reflections.select(nonzero_partiality_sel)
    inverse_partiality.set_selected(
        nonzero_partiality_sel.iselection(), 1.0 / good_refl["partiality"]
    )
    return inverse_partiality


def _calculate_combined_raw_intensities(Iprf, Isum, Vprf, Vsum, Imid):
    """Use partiality-corrected Isum, alongside Iprf to calculate
    combined raw intensities."""
    w = 1.0 / (1.0 + (Isum / Imid) ** 3)
    w.set_selected(Isum <= 0, 1.0)
    Icomb = (w * Iprf) + ((1.0 - w) * Isum)
    Vcomb = (w * Vprf) + ((1.0 - w) * Vsum)
    return Icomb, Vcomb
