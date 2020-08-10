# -*- coding: utf8 -*-

from __future__ import absolute_import, division, print_function

import logging

from cctbx import miller
from cctbx import crystal
from dials.array_family import flex

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")


def compute_mean_weighted_cc_half(intensities):
    cc_bins = intensities.cc_one_half_sigma_tau(use_binning=True, return_n_refl=True)
    bin_data = [b for b in cc_bins.data if b is not None]
    return flex.mean_weighted(
        flex.double(b[0] for b in bin_data), flex.double(b[1] for b in bin_data),
    )


class PerGroupCChalfStatistics(object):
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

        required = ["miller_index", "intensity", "variance", "dataset", "group"]

        for r in required:
            if r not in reflection_table:
                raise KeyError("Column %s not present in reflection table" % r)

        self._intensities = (
            miller.array(
                miller.set(
                    crystal.symmetry(
                        space_group=space_group, unit_cell=mean_unit_cell,
                    ),
                    reflection_table["miller_index"],
                ),
                data=reflection_table["intensity"],
                sigmas=reflection_table["variance"],
            )
            .set_observation_type_xray_intensity()
            .as_non_anomalous_array()
            .map_to_asu()
        )
        sel = self._intensities.resolution_filter_selection(d_min=d_min, d_max=d_max)
        self._intensities = self._intensities.select(sel)
        self._groups = reflection_table["group"].select(sel)
        self._intensities.setup_binner_counting_sorted(n_bins=n_bins)
        self.mean_cchalf = compute_mean_weighted_cc_half(self._intensities)
        logger.info(f"CC 1/2 mean: {self.mean_cchalf:.3f}")
        self.cchalf_i = self._compute_cchalf_excluding_each_group()

    def _compute_cchalf_excluding_each_group(self):
        """Compute the CC½ with each group excluded in turn

        For each group, compute the CC½ excluding reflections in that group.

        Returns (flex.double): The list of CC½ values excluding each group.
        """

        cchalf_i = flex.double()
        for i_group in self._groups.counts():
            intensities = self._intensities.select(self._groups != i_group)
            intensities.use_binning_of(self._intensities)
            cchalf_i.append(compute_mean_weighted_cc_half(intensities))
            logger.info(f"CC½ excluding group {i_group}: {cchalf_i[-1]:.3f}")
        return cchalf_i

    @property
    def delta_cchalf_i(self):
        """Return the ΔCC½ for each group excluded

        Returns (flex.double): The list of ΔCC½ values excluding each group.
        """
        return self.mean_cchalf - self.cchalf_i

    @property
    def group_ids(self):
        """The group ids corresponding to the ΔCC½ values"""
        return flex.size_t(self._groups.counts().keys())
