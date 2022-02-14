from __future__ import annotations

from scitbx.array_family import flex

import dials.extensions
from dials.algorithms.image.filter import convolve


class RadialProfileSpotFinderThresholdExt:
    """
    Extension to calculate a radial profile threshold. This method calculates
    background value and sigma in 2θ shells, then sets a threshold at a level
    n_sigma above the radial background. As such, it is important to have the
    beam centre correct and to mask out any significant shadows. The method may
    be particularly useful for electron diffraction images, where there can be
    considerable inelastic scatter around low resolution spots. In addition, the
    algorithm is relatively insensitive to noise properties of the detector.
    This helps for the case of integrating detectors with poorly known gain
    and response statistics.

    A similar algorithm is available in other programs. The description of
    'peakfinder 8' in https://doi.org/10.1107/S1600576714007626 was helpful
    in the development of this method.
    """

    name = "radial_profile"

    @staticmethod
    def phil():
        from libtbx.phil import parse

        phil = parse(
            """
        n_sigma = 8
          .type = int
          .help = "Sigma multiplier for determining the threshold value"

        blur = narrow wide
          .type = choice
          .help = "Optional preprocessing of the image by a convolution with"
                  "a simple Gaussian kernel of size either 3×3 (narrow) or"
                  "5×5 (wide). This may help to reduce noise peaks and to"
                  "combine split spots."

        n_bins = 100
          .type = int
          .help = "Number of 2θ bins in which to calculate background"
        """
        )
        return phil

    def __init__(self, params):
        """
        Initialise the algorithm.

        :param params: The input parameters
        """
        self.params = params

        # Set approximate Gaussian kernel for blurring
        if self.params.spotfinder.threshold.radial_profile.blur == "narrow":
            # fmt: off
            self.kernel = flex.double(
                (0.0625, 0.125, 0.0625,
                 0.125,  0.25,  0.125,
                 0.0625, 0.125, 0.0625)
            )
            # fmt: on
            self.kernel.reshape(flex.grid((3, 3)))
        elif self.params.spotfinder.threshold.radial_profile.blur == "wide":
            # fmt: off
            self.kernel = (
                flex.double(
                    (
                        1,  4,  7,  4,  1,
                        4, 16, 26, 16,  4,
                        7, 26, 41, 26,  7,
                        4, 16, 26, 16,  4,
                        1,  4,  7,  4,  1,
                    )
                ) / 273
            )
            # fmt: on
            self.kernel.reshape(flex.grid((5, 5)))
        else:
            self.kernel = None

    def compute_threshold(
        self, image, mask, *, imageset, i_panel, region_of_interest=None, **kwargs
    ):
        """
        Compute the threshold.

        :param image: The image to process
        :param mask: The pixel mask on the image
        :**kwargs: Arbitrary keyword arguments
        :returns: A boolean mask showing foreground/background pixels
        """

        if self.kernel:
            image = convolve(image, self.kernel)

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
        slot_info = list(h0.slot_infos())
        lower_bound = slot_info[0].low_cutoff
        lookup = full_two_theta_array - lower_bound
        lookup.set_selected(lookup < 0, 0)

        # Truncate just under the shifted upper bound and rescale
        upper_bound = slot_info[-1].high_cutoff - lower_bound
        lookup.set_selected(lookup >= upper_bound, upper_bound - 1e-10)
        lookup /= upper_bound  # values now in range [0,1)

        # Convert to a size_t lookup into the threshold array
        lookup *= n_bins
        lookup = flex.floor(lookup).iround().as_size_t()

        # Now construct a threshold image
        thresh_im = threshold.select(lookup.as_1d())

        # Peaks are unmasked pixels greater than the threshold
        peaks = image > thresh_im
        return peaks & mask
