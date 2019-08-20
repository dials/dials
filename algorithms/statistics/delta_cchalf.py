#!/usr/bin/env python
#
# delta_cchalf.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function

import logging
from collections import defaultdict
from math import sqrt, floor

from cctbx import miller
from cctbx import crystal, uctbx
from dials.array_family import flex
import six

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")


class ResolutionBinner(object):
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
        dmin_inv_sq = 1.0 / dmin ** 2
        dmax_inv_sq = 1.0 / dmax ** 2
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
                logger.info("%d: %.3f, %.3f" % (i, sqrt(1 / b0), sqrt(1 / b1)))
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
        d2 = 1 / d ** 2
        bin_index = int(floor((d2 - self._xmin) / self._bin_size))
        if bin_index >= self._nbins:
            bin_index = self._nbins - 1
        if bin_index < 0:
            bin_index = 0
        return bin_index


class ReflectionSum(object):
    """
    A helper class to store sums of X and X**2

    """

    def __init__(self, sum_x=0, sum_x2=0, n=0):
        self.sum_x = sum_x
        self.sum_x2 = sum_x2
        self.n = n


class BinData(object):
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
    sigma_y = sum([(m - mean_of_means) ** 2 for m in mean]) / (n - 1)
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


class PerImageCChalfStatistics(object):
    """
    A class to compute per image CC 1/2 statistics

    """

    def __init__(
        self,
        miller_index,
        identifiers,
        dataset,
        images,
        intensity,
        variance,
        unit_cell,
        space_group,
        nbins=10,
        dmin=None,
        dmax=None,
        mode="dataset",
        image_group=10,
    ):
        """
        Initialise

        :param miller_index: The list of miller indices
        :param dataset: The list of dataset numbers
        :param intensity: The list of intensities
        :param variance: The list of variances
        :param unit_cell: The unit cell or list of unit cells
        :param space_group: The space group
        :param nbins: The number of bins
        :param dmin: The maximum resolution
        :param dmax: The minimum resolution

        """

        assert len(set(dataset)) == len(unit_cell)

        # Reject reflections with negative variance
        selection = variance > 0
        miller_index = miller_index.select(selection)
        dataset = dataset.select(selection)
        intensity = intensity.select(selection)
        variance = variance.select(selection)
        images = images.select(selection)

        # Compute mean unit_cell
        if len(unit_cell) == 1:
            mean_unit_cell = unit_cell[0]
        else:
            mean_parameters = [0, 0, 0, 0, 0, 0]
            for uc in unit_cell:
                for i in range(6):
                    mean_parameters[i] += uc.parameters()[i]
            for i in range(6):
                mean_parameters[i] /= len(unit_cell)
            mean_unit_cell = uctbx.unit_cell(mean_parameters)

        # Map the miller indices to the ASU
        D = flex.double(len(dataset), 0)
        for i, uc in zip(identifiers, unit_cell):

            # Select reflections
            selection = dataset == i
            hkl = miller_index.select(selection)

            # Compute resolution
            d = flex.double([uc.d(h) for h in hkl])
            D.set_selected(selection, d)

            # Compute asu miller index
            cs = crystal.symmetry(uc, space_group=space_group)
            ms = miller.set(cs, hkl)
            ms_asu = ms.map_to_asu()
            miller_index.set_selected(selection, ms_asu.indices())

        assert all(d > 0 for d in D)

        # Filter by dmin and dmax
        if dmin is not None:
            selection = D > dmin
            miller_index = miller_index.select(selection)
            dataset = dataset.select(selection)
            intensity = intensity.select(selection)
            variance = variance.select(selection)
            D = D.select(selection)
            images = images.select(selection)

        if dmax is not None:
            selection = D < dmax
            miller_index = miller_index.select(selection)
            dataset = dataset.select(selection)
            intensity = intensity.select(selection)
            variance = variance.select(selection)
            D = D.select(selection)
            images = images.select(selection)

        # Save the arrays
        self._miller_index = miller_index
        self._dataset = dataset
        self._intensity = intensity
        self._variance = variance

        # Create the resolution bins
        self._num_bins = nbins
        self._dmin = dmin
        self._dmax = dmax
        if dmin is None:
            self._dmin = min(D)
        if dmax is None:
            self._dmax = max(D)
        binner = ResolutionBinner(mean_unit_cell, self._dmin, self._dmax, nbins)

        # Create lookups for elements by miller index
        index_lookup = defaultdict(list)
        for i, h in enumerate(miller_index):
            index_lookup[h].append(i)

        # Compute the Overall Sum(X) and Sum(X^2) for each unique reflection
        reflection_sums = defaultdict(ReflectionSum)
        for h in index_lookup:
            I = [intensity[i] for i in index_lookup[h]]
            n = len(I)
            sum_x = sum(I)
            sum_x2 = sum(II ** 2 for II in I)
            reflection_sums[h] = ReflectionSum(sum_x, sum_x2, n)

        # Compute some numbers
        self._num_datasets = len(set(dataset))
        self._num_reflections = len(miller_index)
        self._num_unique = len(reflection_sums)

        logger.info("")
        logger.info("# Datasets: %s" % self._num_datasets)
        logger.info("# Reflections: %s" % self._num_reflections)
        logger.info("# Unique: %s" % self._num_unique)

        # Compute the CC 1/2 for all the data
        self._cchalf_mean = self._compute_cchalf(reflection_sums, binner)
        logger.info("CC 1/2 mean: %.3f" % (100 * self._cchalf_mean))

        # override dataset here with a batched-dependent
        self.expid_to_image_groups = {id_: [] for id_ in set(dataset)}
        if mode == "image_group":
            image_groups = flex.int(dataset.size(), 0)
            self.image_group_to_expid_and_range = {}

            counter = 0
            for id_ in set(dataset):
                sel = dataset == id_
                images_in_dataset = images.select(sel)
                unique_images = set(images_in_dataset)
                min_img, max_img = (min(unique_images), max(unique_images))
                min_img = (int(floor(min_img / image_group)) * image_group) + 1
                for i in range(min_img, max_img + 1, image_group):
                    group_sel = (images_in_dataset >= i) & (
                        images_in_dataset < i + image_group
                    )
                    image_groups.set_selected(
                        (sel.iselection().select(group_sel)), counter
                    )
                    self.image_group_to_expid_and_range[counter] = (
                        id_,
                        (i, i + image_group - 1),
                    )
                    self.expid_to_image_groups[id_].append(counter)
                    counter += 1
            self._cchalf = self._compute_cchalf_excluding_each_dataset(
                reflection_sums, binner, miller_index, image_groups, intensity
            )

        else:
            self._cchalf = self._compute_cchalf_excluding_each_dataset(
                reflection_sums, binner, miller_index, dataset, intensity
            )

    def _compute_cchalf(self, reflection_sums, binner):
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

    def _compute_cchalf_excluding_each_dataset(
        self, reflection_sums, binner, miller_index, dataset, intensity
    ):
        """
        Compute the CC 1/2 with an image excluded.

        For each image, update the sums by removing the contribution from the image
        and then compute the CC 1/2 of the remaining data

        """

        # Create a lookup table for each reflection by dataset
        dataset_lookup = defaultdict(list)
        for i, b in enumerate(dataset):
            dataset_lookup[b].append(i)

        # Compute CC1/2 minus each dataset
        cchalf_i = {}
        for dataset in dataset_lookup:

            # Find all observations from this dataset and create a lookup based on
            # miller index
            index_lookup = defaultdict(list)
            for i in dataset_lookup[dataset]:
                index_lookup[miller_index[i]].append(i)

            # Loop through all the reflections and remove the contribution of
            # reflections from the current dataset
            dataset_reflection_sums = defaultdict(ReflectionSum)
            for h in reflection_sums:
                sum_x = reflection_sums[h].sum_x
                sum_x2 = reflection_sums[h].sum_x2
                n = reflection_sums[h].n
                for i in index_lookup[h]:
                    sum_x -= self._intensity[i]
                    sum_x2 -= self._intensity[i] ** 2
                    n -= 1
                dataset_reflection_sums[h] = ReflectionSum(sum_x, sum_x2, n)

            # Compute the CC 1/2 without the reflections from the current dataset
            cchalf = self._compute_cchalf(dataset_reflection_sums, binner)
            cchalf_i[dataset] = cchalf
            logger.info("CC 1/2 excluding dataset %d: %.3f" % (dataset, 100 * cchalf))

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
        Return the Delta CC 1/2 for each image excluded

        """
        return {k: self._cchalf_mean - v for k, v in six.iteritems(self._cchalf)}
