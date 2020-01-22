"""
An augmentation of iotbx.merging_statistics, kinder in the case of few reflections.

This module provides a decorated version of iotbx.merging_statistics.data_statistics
that tries to be more accommodating in the case of few reflections, as in a small
molecule or small wedge data set.
"""

from __future__ import absolute_import, division, print_function

from iotbx import merging_statistics


def dataset_statistics(*args, **kwargs):
    # type: (...) -> merging_statistics.dataset_statistics
    """
    Calculate merging statistics, handling errors arising from too few reflections.

    Try to return a iotbx.merging_statistics.dataset_statistics instance with default
    parameters.  If a StatisticsError is encountered because of too few reflections,
    try again using the 'counting_sorted' binning algorithm.  If that still doesn't
    work, progressively reduce the number of bins until successful, or re-raise the
    StatisticsError if no success is possible with fewer than six bins.

    Args:
        *args:  Positional arguments of iotbx.merging_statistics.dataset_statistics.
        **kwargs:  Keyword arguments of iotbx.merging_statistics.dataset_statistics.

    Returns:
        An iotbx. merging_statistics.dataset_statistics object containing the
        calculated merging statistics.
    """
    while True:
        try:
            return merging_statistics.dataset_statistics(*args, **kwargs)
        except merging_statistics.StatisticsError:
            if kwargs["binning_method"] != "counting_sorted":
                kwargs["binning_method"] = "counting_sorted"
            elif kwargs["n_bins"] > 8:
                kwargs["n_bins"] -= 3
            else:
                raise

            continue
