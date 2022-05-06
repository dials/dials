from __future__ import annotations

import logging

import libtbx

from dials.algorithms.spot_finding.threshold import (
    DispersionZeroTruncatedThresholdStrategy,
)
from dials.extensions.dispersion_extended_spotfinder_threshold_ext import (
    estimate_global_threshold,
)

logger = logging.getLogger(__name__)


class DispersionZeroTruncatedSpotFinderThresholdExt:
    """Extensions to do dispersion threshold."""

    name = "dispersion_zero_truncated"

    default = False

    @staticmethod
    def phil():
        return None

    def __init__(self, params):
        """
        Initialise the algorithm.

        :param params: The input parameters
        """
        self.params = params

    def compute_threshold(self, image, mask, **kwargs):
        r"""
        Compute the threshold.

        :param image: The image to process
        :param mask: The pixel mask on the image
        :\*\*kwargs: Arbitrary keyword arguments
        :returns: A boolean mask showing foreground/background pixels
        """

        params = self.params
        if params.spotfinder.threshold.dispersion.global_threshold is libtbx.Auto:
            params.spotfinder.threshold.dispersion.global_threshold = int(
                estimate_global_threshold(image, mask)
            )
            logger.info(
                "Setting global_threshold: %i",
                params.spotfinder.threshold.dispersion.global_threshold,
            )

        self._algorithm = DispersionZeroTruncatedThresholdStrategy(
            kernel_size=params.spotfinder.threshold.dispersion.kernel_size,
            gain=params.spotfinder.threshold.dispersion.gain,
            mask=params.spotfinder.lookup.mask,
            n_sigma_b=params.spotfinder.threshold.dispersion.sigma_background,
            n_sigma_s=params.spotfinder.threshold.dispersion.sigma_strong,
            min_count=params.spotfinder.threshold.dispersion.min_local,
            global_threshold=params.spotfinder.threshold.dispersion.global_threshold,
        )

        return self._algorithm(image, mask)
