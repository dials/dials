#!/usr/bin/env python
#
#  mcd.py
#
#  Copyright (C) 2015 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from dials.algorithms.refinement.outlier_detection import CentroidOutlier
from dials.algorithms.statistics.fast_mcd import FastMCD, maha_dist_sq
from dials_refinement_helpers_ext import qchisq


class MCD(CentroidOutlier):
    """Implementation of the CentroidOutlier class using Mahalanobis distance
    calculated using the robust location and scatter estimate from the Minimum
    Covariance Determinant estimate."""

    def __init__(
        self,
        cols=None,
        min_num_obs=20,
        separate_experiments=True,
        separate_panels=True,
        block_width=None,
        alpha=0.5,
        max_n_groups=5,
        min_group_size=300,
        n_trials=500,
        k1=2,
        k2=2,
        k3=100,
        threshold_probability=0.975,
    ):

        if cols is None:
            cols = ["x_resid", "y_resid", "phi_resid"]
        CentroidOutlier.__init__(
            self,
            cols=cols,
            min_num_obs=min_num_obs,
            separate_experiments=separate_experiments,
            separate_panels=separate_panels,
            block_width=block_width,
        )

        # Keep the FastMCD options here
        self._alpha = alpha
        self._max_n_groups = max_n_groups
        self._min_group_size = min_group_size
        self._n_trials = n_trials
        self._k1 = k1
        self._k2 = k2
        self._k3 = k3

        # Calculate Mahalanobis distance threshold
        df = len(cols)
        self._mahasq_cutoff = qchisq(threshold_probability, df)

        return

    def _detect_outliers(self, cols):

        fast_mcd = FastMCD(
            cols,
            alpha=self._alpha,
            max_n_groups=self._max_n_groups,
            min_group_size=self._min_group_size,
            n_trials=self._n_trials,
            k1=self._k1,
            k2=self._k2,
            k3=self._k3,
        )

        # get location and MCD scatter estimate
        T, S = fast_mcd.get_corrected_T_and_S()

        # get squared Mahalanobis distances
        d2s = maha_dist_sq(cols, T, S)

        # compare to the threshold
        outliers = d2s > self._mahasq_cutoff

        return outliers
