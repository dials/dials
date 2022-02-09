from __future__ import annotations

from scitbx.array_family import flex

import dials.extensions


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

        n_sigma = 6
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

        # Initial dispersion threshold to determine likely background pixels
        Dispersion = dials.extensions.SpotFinderThreshold.load("dispersion")
        threshold_function = Dispersion(self.params)
        peak_pixels = threshold_function.compute_threshold(image, mask)
        background_pixels = mask & ~peak_pixels
        background = image.select(background_pixels.iselection())

        # compute histogram of two-theta values, then same weighted
        # by pixel values, finally divide latter by former to get
        # the radial profile out, need to set the number of bins
        # sensibly; inspired by method in PyFAI
        panel = imageset.get_detector()[i_panel]
        beam = imageset.get_beam()

        # Get 2θ array for the background pixels
        full_two_theta_array = panel.get_two_theta_array(beam.get_s0())
        if region_of_interest:
            x0, x1, y0, y1 = region_of_interest
            full_two_theta_array = full_two_theta_array[y0:y1, x0:x1]
        two_theta_array = full_two_theta_array.as_1d().select(
            background_pixels.iselection()
        )

        # Use flex.weighted_histogram
        n_bins = self.params.spotfinder.threshold.radial_profile.n_bins
        h0 = flex.weighted_histogram(two_theta_array, n_slots=n_bins)
        h1 = flex.weighted_histogram(two_theta_array, background, n_slots=n_bins)
        h2 = flex.weighted_histogram(
            two_theta_array, background * background, n_slots=n_bins
        )

        d0 = h0.slots()
        d1 = h1.slots()
        d2 = h2.slots()

        I = d1 / d0
        I2 = d2 / d0
        sig = flex.sqrt(I2 - flex.pow2(I))

        # Determine the threshold value for each bin
        n_sigma = self.params.spotfinder.threshold.radial_profile.n_sigma
        threshold = I + n_sigma * sig

        # Shift the full 2θ array to the lower bound and truncate
        infos = list(h0.slot_infos())
        lower_bound = infos[0].low_cutoff
        lookup = full_two_theta_array - lower_bound
        lookup.set_selected(lookup < 0, 0)

        # Truncate just under the shifted upper bound and rescale
        upper_bound = infos[-1].high_cutoff - lower_bound
        lookup.set_selected(lookup >= upper_bound, upper_bound - 1e-10)
        lookup /= upper_bound  # values now in range [0,1)

        # Convert to a size_t lookup into the threshold array
        lookup *= n_bins
        lookup = flex.floor(lookup).iround().as_size_t()

        # Now construct a threshold image
        thresh_im = threshold.select(lookup.as_1d())

        # Peaks are greater than the threshold
        peaks = image > thresh_im
        return peaks
