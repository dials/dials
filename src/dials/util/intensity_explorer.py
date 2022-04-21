"""
Examine the distribution of diffraction spot intensities.

This module defines a class IntensityDist, with several methods for exploring
the distribution of measured spot intensities in an X-ray diffraction
experiment.  The user may wish to use this information to inform decisions
regarding the error model employed in analysing the data.  Data are passed in
as an unmerged MTZ file (see http://www.ccp4.ac.uk/html/mtzformat.html) and the
resulting IntensityDist instance contains the pertinent columns of data, along
with normal order statistic medians of the z-scores of the intensities, for
constructing a normal probability plot (See
https://www.itl.nist.gov/div898/handbook/eda/section3/normprpl.htm).
"""

from __future__ import annotations

import logging

import scipy.stats

from cctbx import miller
from dxtbx.model import ExperimentList

from dials.array_family import flex

log = logging.getLogger("dials.util.intensity_explorer")


class IntensityDist:
    def __init__(
        self,
        rtable,
        elist,
        calculate_variances=False,
        keep_singles=False,
        uncertainty="sigma",
        outfile=None,
    ):
        """
        Generate z-scores and a normal probability plot from a DIALS
        reflection_table and a dxtbx ExperimentList, containing the observations
        and the corresponding experiments, respectively.

        :param rtable: A reflection table object, containing at least the columns
          * ``miller_index``
          * ``intensity.sum.value``
          * ``intensity.sum.variance``
          * ``xyzobs.px.value``
        :type rtable: dials.array_family_flex_ext.reflection_table
        :param elist: A corresponding experiment list.
        :type elist: dxtbx.model.ExperimentList
        :param calculate_variances: Choose whether to calculate weighted
        aggregate variances.  Doing so incurs a performance penalty.
        Defaults to False.
        :type calculate_variances: bool
        :param keep_singles: Choose whether to keep multiplicity-1 reflections.
        Defaults to False.
        :type keep_singles: bool
        :param uncertainty: Measure of spread to use in normalising the
        z-scores, i.e. z = (I - <I>) / uncertainty.
        Possible values for uncertainty:
        * 'sigma':    Use measured sigma values;
        * 'stddev':   Use sample standard deviations calculated as
                      square-root of unbiased weighted sample variances
                      of symmetry-equivalent reflection intensities;
        Defaults to 'sigma'.
        :type uncertainty: str
        :param outfile: Filename root for output PNG plots.
        Defaults to None.
        :type: outfile: str
        """

        if not isinstance(rtable, flex.reflection_table) or not isinstance(
            elist, ExperimentList
        ):
            raise TypeError(
                "Must be called with a reflection table and an experiment list."
            )

        rtable = rtable.copy()
        # Discard unindexed reflections (only necessary because of
        # https://github.com/dials/dials/issues/615 —
        # TODO remove the line below when issue #615 is fixed).
        rtable.del_selected(rtable["id"] == -1)
        rtable["miller_index.asu"] = rtable["miller_index"]

        # Divide reflections by the space group to which they have been indexed.
        self.rtables = {
            expt.crystal.get_space_group().make_tidy(): flex.reflection_table()
            for expt in elist
        }
        for expt, sel in rtable.iterate_experiments_and_indices(elist):
            sg = expt.crystal.get_space_group().make_tidy()
            self.rtables[sg].extend(rtable.select(sel))
        # Map Miller indices to asymmetric unit.
        for space_group, rtable in self.rtables.items():
            # TODO Handle anomalous flag sensibly.  Currently assumes not anomalous.
            miller.map_to_asu(space_group.type(), False, rtable["miller_index.asu"])

        # Calculate normal probability plot data.
        self._multiplicity_mean_error_stddev(
            calculate_variances=calculate_variances, keep_singles=keep_singles
        )
        self._make_z(uncertainty)
        self._probplot_data()

        self.rtable = flex.reflection_table()
        for rtable in self.rtables.values():
            self.rtable.extend(rtable)

        if not outfile:
            outfile = ""
        self.outfile = outfile

    def _multiplicity_mean_error_stddev(
        self, calculate_variances=False, keep_singles=False
    ):
        """
        Calculate aggregate properties of grouped symmetry-equivalent reflections.

        Populate the reflection table of observations with the following
        properties:
          * ``multiplicity`` — Multiplicity of observations of a given reflection
          in the asymmetric unit;
          :type: `dials.array_family_flex_ext.int` array
          * ``intensity.mean.value`` — Mean of symmetry-equivalent reflections,
          weighted by measurement error;
          :type: `dials.array_family_flex_ext.double` array
          * ``intensity.mean.std_error`` — Standard error on the weighted mean;
          :type: `dials.array_family_flex_ext.double` array
          * (optional) ``intensity.mean.variance`` — variance of
          symmetry-equivalent reflections, weighted by measurement error;
          :type: `dials.array_family_flex_ext.double` array

        :param calculate_variances: Elect whether to calculate the weighted
        variances.  Defaults to False, to spare an expensive computation.
        :type calculate_variances: bool
        :param keep_singles: Choose whether to keep single-multiplicity
        reflections.
        :type keep_singles: bool
        """

        for key, rtable in self.rtables.items():
            # Sort the reflection table for speedier iteration.
            rtable.sort("miller_index.asu")
            # Record the positions of any multiplicity-1 reflections.
            if not keep_singles:
                singles = flex.size_t()
            # Record the multiplicities.
            multiplicity = flex.int()
            # For weighted averaging.
            weights = 1 / rtable["intensity.sum.variance"]
            sum_weights = flex.double()
            if calculate_variances:
                sum_square_weights = flex.double()
            # Calculate the weighted mean intensities.
            i_means = flex.double()
            # Calculate the standard deviations from unbiased weighted variances.
            variances = flex.double()

            # Iterate over the reflections, grouping by equivalent Miller index,
            # to calculate multiplicities, weighted mean intensities, etc..
            # Some time can be saved by only calculating variances if necessary.
            # Initial values:
            prev_index = None
            count = 1
            # The following will be set during loop iteration
            i_sum, sum_weight, sum_square_weight = None, None, None
            # One big loop through the entire reflection table:
            for j in range(rtable.size()):
                index = rtable["miller_index.asu"][j]
                weight = weights[j]
                # Aggregate within a symmetry-equivalent group of reflections:
                if index == prev_index:
                    count += 1
                    i_sum += weight * rtable["intensity.sum.value"][j]
                    sum_weight += weight
                    if calculate_variances:
                        sum_square_weight += weight * weight
                # Record the aggregated values for the group:
                elif prev_index:
                    if count == 1 and not keep_singles:
                        singles.append(j - 1)
                    multiplicity.extend(flex.int(count, count))
                    i_means.extend(flex.double(count, i_sum / sum_weight))
                    sum_weights.extend(flex.double(count, sum_weight))
                    if calculate_variances:
                        sum_square_weights.extend(flex.double(count, sum_square_weight))
                    # And reinitialise:
                    prev_index = index
                    count = 1
                    i_sum = weight * rtable["intensity.sum.value"][j]
                    sum_weight = weight
                    if calculate_variances:
                        sum_square_weight = weight * weight
                # Handle the first row:
                else:
                    prev_index = rtable["miller_index.asu"][j]
                    i_sum = weight * rtable["intensity.sum.value"][j]
                    sum_weight = weight
                    if calculate_variances:
                        sum_square_weight = weight * weight
            # Record the aggregated values for the last group:
            if count == 1 and not keep_singles:
                singles.append(rtable.size() - 1)
            multiplicity.extend(flex.int(count, count))
            i_means.extend(flex.double(count, i_sum / sum_weight))
            sum_weights.extend(flex.double(count, sum_weight))
            if calculate_variances:
                sum_square_weights.extend(flex.double(count, sum_square_weight))

            # Discard singletons:
            if not keep_singles:
                singles_del = flex.bool(rtable.size(), True)
                singles_del.set_selected(singles, False)
                multiplicity, weights, sum_weights, i_means = [
                    a.select(singles_del)
                    for a in (multiplicity, weights, sum_weights, i_means)
                ]
                rtable.del_selected(singles)
                if calculate_variances:
                    sum_square_weights = sum_square_weights.select(singles_del)

            # Record the multiplicities in the reflection table.
            rtable["multiplicity"] = multiplicity
            # Record the weighted mean intensities in the reflection table.
            rtable["intensity.mean.value"] = i_means
            # Record the standard errors on the means in the reflection table.
            rtable["intensity.mean.std_error"] = flex.sqrt(1 / sum_weights)

            if calculate_variances:
                # Initialise values:
                prev_index = None
                weighted_sum_square_residual = None
                for j in range(rtable.size()):
                    index = rtable["miller_index.asu"][j]
                    weight = weights[j]
                    residual = rtable["intensity.sum.value"][j] - i_means[j]
                    # Aggregate within a symmetry-equivalent group of reflections:
                    if index == prev_index:
                        count += 1
                        weighted_sum_square_residual += weight * residual * residual
                    # Record the aggregated value for the group:
                    elif prev_index:
                        # The weighted variance is undefined for multiplicity=1,
                        # use the measured variance instead in this case.
                        if count == 1:
                            variances.append(rtable["intensity.sum.variance"][j - 1])
                        else:
                            sum_weight = sum_weights[j - 1]
                            var_weight = 1 / (
                                sum_weight - sum_square_weights[j - 1] / sum_weight
                            )
                            variances.extend(
                                flex.double(
                                    count, weighted_sum_square_residual * var_weight
                                )
                            )
                        # Reinitialise:
                        prev_index = index
                        count = 1
                        weighted_sum_square_residual = weight * residual * residual
                    # Handle the first row:
                    else:
                        prev_index = rtable["miller_index.asu"][j]
                        count = 1
                        weighted_sum_square_residual = weight * residual * residual
                # Record the aggregated values for the last group:
                # The weighted variance is undefined for multiplicity=1,
                # use the measured variance instead in this case.
                if count == 1:
                    variances.append(rtable["intensity.sum.variance"][-1])
                else:
                    sum_weight = sum_weights[-1]
                    var_weight = 1 / (sum_weight - sum_square_weights[-1] / sum_weight)
                    variances.extend(
                        flex.double(count, weighted_sum_square_residual * var_weight)
                    )
                # Record the variances in the reflection table.
                rtable["intensity.mean.variance"] = variances

            self.rtables[key] = rtable

    def _make_z(self, uncertainty="sigma"):
        """
        Generate reflection z-scores.

        Calculate z-scores from reflection intensities, weighted mean
        intensities and a chosen measure of uncertainty in the intensity
        measurement.

        :param uncertainty: Chosen measure of uncertainty.  Options are
          * ``stddev`` — standard deviation, as calculated from the unbiased
          weighted variance aggregated amongst all symmetry-equivalent reflections;
          * ``sigma`` — measurement error for individual reflections.
        :type uncertainty: str
        """

        uncertainty_name = {
            "sigma": "intensity.sum.variance",
            "stddev": "intensity.mean.variance",
        }[uncertainty]

        for key, rtable in self.rtables.items():
            try:
                uncertainty_value = flex.sqrt(rtable[uncertainty_name])
            except KeyError:
                uncertainty_value = flex.sqrt(rtable["intensity.sum.variance"])
                log.warn(
                    """Weighted variances haven't been calculated,
      be sure to specify calculate_variances=True to use them.
      Defaulting to measured σ values as a measure of uncertainty instead."""
                )

            z = (
                rtable["intensity.sum.value"] - rtable["intensity.mean.value"]
            ) / uncertainty_value
            rtable["intensity.z_score"] = z

            self.rtables[key] = rtable

    def _probplot_data(self):
        """Generate the data for a normal probability plot of z-scores."""

        for key, rtable in self.rtables.items():
            order = flex.sort_permutation(rtable["intensity.z_score"])
            osm = flex.double(rtable.size(), 0)
            probplot = scipy.stats.probplot(rtable["intensity.z_score"], fit=False)
            osm.set_selected(order, flex.double(probplot[0]))
            rtable["intensity.order_statistic_medians"] = osm

            self.rtables[key] = rtable
