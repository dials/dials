from __future__ import annotations

import logging

from scitbx.array_family import flex

logger = logging.getLogger("dials.extensions.radial_profile_spotfinder_threshold_ext")


class RadialProfileSpotFinderThresholdExt:
    """Extension to do radial profile thresholding."""

    name = "radial_profile"

    @staticmethod
    def phil():
        from libtbx.phil import parse

        phil = parse(
            """
        n_bins = 100
        .type = int

        n_sigma = 3
        .type = int
        .help = "Sigma multiplier for determining the threshold value"
    """
        )
        return phil

    def __init__(self, params):
        """
        Initialise the algorithm.

        :param params: The input parameters
        """
        self.params = params

    def compute_threshold(self, image, mask, **kwargs):
        """
        Compute the threshold.

        :param image: The image to process
        :param mask: The pixel mask on the image
        :**kwargs: Arbitrary keyword arguments
        :returns: A boolean mask showing foreground/background pixels
        """

        try:
            imageset = kwargs["imageset"]
            i_panel = kwargs["i_panel"]
        except KeyError:
            raise ValueError("Missing required arguments")
        region_of_interest = kwargs.get("region_of_interest")

        # TODO - algorithm goes here
        result = flex.bool(len(image))
        result.reshape(image.accessor())

        return result
