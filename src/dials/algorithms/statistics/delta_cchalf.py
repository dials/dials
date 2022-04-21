from __future__ import annotations

import logging
from collections import defaultdict
from math import floor, sqrt

from cctbx import crystal, miller

from dials.array_family import flex

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")


class ResolutionBinner:
    """
    A class to bin the data by resolution
    """

    def __init__(self, unit_cell, dmin, dmax, nbins, output=True):
        """
        Initialise the binner

        :param unit_cell: A unit cell object
        :param dmin: The maximum resolution
        :param dmax: The minimum resolution
        :param nbins: The number of bins
        """
        if output:
            logger.info("Resolution bins")
        assert dmin < dmax
        dmin_inv_sq = 1.0 / dmin**2
        dmax_inv_sq = 1.0 / dmax**2
        self._unit_cell = unit_cell
        self._nbins = nbins
        self._xmin = dmax_inv_sq
        self._xmax = dmin_inv_sq
        assert self._xmin < self._xmax
        self._bin_range = self._xmax - self._xmin
        self._bin_size = self._bin_range / self._nbins
        self._bins = []
        for i in range(self._nbins):
            b0, b1 = (
                self._xmin + i * self._bin_size,
                self._xmin + (i + 1) * self._bin_size,
            )
            if output:
                logger.info("%d: %.3f, %.3f", i, sqrt(1 / b0), sqrt(1 / b1))
            self._bins.append((b0, b1))

    def nbins(self):
        """
        :returns: The number of bins
        """
        return self._nbins

    def index(self, h):
        """
        Get the bin index from the miller index

        :param h: The miller index
        :returns: The bin index
        """
        d = self._unit_cell.d(h)
        d2 = 1 / d**2
        bin_index = int(floor((d2 - self._xmin) / self._bin_size))
        if bin_index >= self._nbins:
            bin_index = self._nbins - 1
        if bin_index < 0:
            bin_index = 0
        return bin_index


class ReflectionSum:
    """
    A helper class to store sums of X and X**2
    """

    def __init__(self, sum_x=0, sum_x2=0, n=0):
        self.sum_x = sum_x
        self.sum_x2 = sum_x2
        self.n = n


class BinData:
    """
    A helper class to store mean and variance
    """

    def __init__(self):
        self.mean = []
        self.var = []


def compute_cchalf(mean, var):
    """
    Compute the CC 1/2 using the formular from Assmann, Brehm and Diederichs 2016

    :param mean: The list of mean intensities
    :param var: The list of variances on the half set of mean intensities
    :returns: The CC 1/2
    """
    assert len(mean) == len(var)
    n = len(mean)
    mean_of_means = sum(mean) / n
    sigma_e = sum(var) / n
    sigma_y = sum((m - mean_of_means) ** 2 for m in mean) / (n - 1)
    cchalf = (sigma_y - sigma_e) / (sigma_y + sigma_e)
    return cchalf


def compute_mean_cchalf_in_bins(bin_data):
    """
    Compute the mean cchalf averaged across resolution bins

    :param bin_data: The mean and variance in each bin
    :returns: The mean CC 1/2
    """
    mean_cchalf = 0
    count = 0
    for bin_i in bin_data:
        mean = bin_i.mean
        var = bin_i.var
        n = len(mean)
        if n > 1:
            cchalf = compute_cchalf(mean, var)
            mean_cchalf += n * cchalf
            count += n
    mean_cchalf /= count
    return mean_cchalf


def compute_cchalf_from_reflection_sums(reflection_sums, binner):
    """
    Compute the CC 1/2 by computing the CC 1/2 in resolution bins and then
    computing the weighted mean of the binned CC 1/2 values
    """
    # Compute Mean and variance of reflection intensities
    bin_data = [BinData() for _ in range(binner.nbins())]
    for h in reflection_sums:
        sum_x = reflection_sums[h].sum_x
        sum_x2 = reflection_sums[h].sum_x2
        n = reflection_sums[h].n

        if n > 1:
            mean = sum_x / n
            var = (sum_x2 - (sum_x) ** 2 / n) / (n - 1)
            var = var / n
            index = binner.index(h)
            bin_data[index].mean.append(mean)
            bin_data[index].var.append(var)

    # Compute the mean cchalf in resolution bins
    return compute_mean_cchalf_in_bins(bin_data)


class PerGroupCChalfStatistics:
    def __init__(
        self,
        reflection_table,
        mean_unit_cell,
        space_group,
        d_min=None,
        d_max=None,
        n_bins=10,
    ):
        # here dataset is the sweep number, group is the group number for doing
        # the cc half analysis. May be the same as dataset if doing per dataset
        # stats, but can also be different if doing image groups.

        self.mean_unit_cell = mean_unit_cell
        self.space_group = space_group
        self.d_min = d_min
        self.d_max = d_max
        self._num_bins = n_bins
        self.reflection_table = reflection_table
        self._cchalf_mean = None
        self._cchalf = None

        required = ["miller_index", "intensity", "variance", "dataset", "group"]

        for r in required:
            if r not in reflection_table:
                raise KeyError(f"Column {r} not present in reflection table")
        # now do a prefiltering
        self.map_to_asu()
        self.d_filter()

        # Create the resolution bins
        if self.d_min is None:
            self.d_min = flex.min(self.reflection_table["d"])
        if self.d_max is None:
            self.d_max = flex.max(self.reflection_table["d"])
        self.binner = ResolutionBinner(mean_unit_cell, self.d_min, self.d_max, n_bins)

        self.reflection_sums = defaultdict(ReflectionSum)
        self.compute_overall_stats()

    def compute_overall_stats(self):
        # Create lookups for elements by miller index
        index_lookup = defaultdict(list)
        for i, h in enumerate(self.reflection_table["miller_index"]):
            index_lookup[h].append(i)

        # Compute the Overall Sum(X) and Sum(X^2) for each unique reflection

        for h in index_lookup:
            sel = flex.size_t(index_lookup[h])
            intensities = self.reflection_table["intensity"].select(sel)
            n = intensities.size()
            sum_x = flex.sum(intensities)
            sum_x2 = flex.sum(flex.pow2(intensities))
            self.reflection_sums[h] = ReflectionSum(sum_x, sum_x2, n)

        # Compute some numbers
        self._num_datasets = len(set(self.reflection_table["dataset"]))
        self._num_groups = len(set(self.reflection_table["group"]))
        self._num_reflections = self.reflection_table.size()
        self._num_unique = len(self.reflection_sums)

        logger.info(
            """
Summary of input data:
# Datasets: %s
# Groups: %s
# Reflections: %s
# Unique reflections: %s""",
            self._num_datasets,
            self._num_groups,
            self._num_reflections,
            self._num_unique,
        )

    def map_to_asu(self):
        """Map the miller indices to the ASU"""
        # d = flex.double([self.mean_unit_cell.d(h) for h in reflection_table["miller_index"]])
        cs = crystal.symmetry(self.mean_unit_cell, space_group=self.space_group)
        ms = miller.set(cs, self.reflection_table["miller_index"], anomalous_flag=False)
        ms_asu = ms.map_to_asu()
        self.reflection_table["miller_index"] = ms_asu.indices()
        self.reflection_table["d"] = ms_asu.d_spacings().data()

    def d_filter(self):
        """Filter on d_min, d_max"""
        if (not self.d_min) and (not self.d_max):
            return
        if self.d_min is not None:
            sel = self.reflection_table["d"] > self.d_min
            if self.d_max is not None:
                sel &= self.reflection_table["d"] < self.d_max
        else:
            sel = self.reflection_table["d"] < self.d_max
        self.reflection_table.select(sel)

    def run(self):
        """Compute the ΔCC½ for all the data"""
        self._cchalf_mean = compute_cchalf_from_reflection_sums(
            self.reflection_sums, self.binner
        )
        logger.info("CC 1/2 mean: %.3f", (100 * self._cchalf_mean))
        self._cchalf = self._compute_cchalf_excluding_each_group()

    def _compute_cchalf_excluding_each_group(self):
        #   self, reflection_sums, binner, miller_index, dataset, intensity
        # ):
        """
        Compute the CC 1/2 with an image excluded.

        For each image, update the sums by removing the contribution from the image
        and then compute the CC 1/2 of the remaining data
        """

        # Create a lookup table for each reflection by dataset i.e. group of images
        dataset_lookup = defaultdict(list)
        for i, b in enumerate(self.reflection_table["group"]):
            dataset_lookup[b].append(i)

        # Compute CC1/2 minus each dataset
        cchalf_i = {}
        for dataset in dataset_lookup:

            # Find all observations from this dataset and create a lookup based on
            # miller index
            index_lookup = defaultdict(list)
            for i in dataset_lookup[dataset]:
                index_lookup[self.reflection_table["miller_index"][i]].append(i)

            # Loop through all the reflections and remove the contribution of
            # reflections from the current dataset
            dataset_reflection_sums = defaultdict(ReflectionSum)
            for h in self.reflection_sums:
                sum_x = self.reflection_sums[h].sum_x
                sum_x2 = self.reflection_sums[h].sum_x2
                n = self.reflection_sums[h].n
                for i in index_lookup[h]:
                    sum_x -= self.reflection_table["intensity"][i]
                    sum_x2 -= self.reflection_table["intensity"][i] ** 2
                    n -= 1
                dataset_reflection_sums[h] = ReflectionSum(sum_x, sum_x2, n)

            # Compute the CC 1/2 without the reflections from the current dataset
            cchalf = compute_cchalf_from_reflection_sums(
                dataset_reflection_sums, self.binner
            )
            cchalf_i[dataset] = cchalf
            logger.info("CC 1/2 excluding group %d: %.3f", dataset, 100 * cchalf)

        return cchalf_i

    def num_datasets(self):
        """
        Return the number of datasets
        """
        return len(self._cchalf)

    def num_reflections(self):
        """
        Return the number of reflections
        """
        return self._num_reflections

    def num_unique(self):
        """
        Return the number of unique reflections
        """
        return self._num_unique

    def mean_cchalf(self):
        """
        Return the mean CC 1/2
        """
        return self._cchalf_mean

    def cchalf_i(self):
        """
        Return the CC 1/2 for each image excluded
        """
        return self._cchalf

    def delta_cchalf_i(self):
        """
        Return the ΔCC½ for each image excluded
        """
        return {k: self._cchalf_mean - v for k, v in self._cchalf.items()}
