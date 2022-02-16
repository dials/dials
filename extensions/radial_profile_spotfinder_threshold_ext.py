from __future__ import annotations

from scitbx.array_family import flex

from dials.algorithms.image.filter import convolve
from dials.algorithms.statistics import BinnedStatistics

# Module-level definition imported by the image viewer
phil_str = """
n_iqr = 6
    .type = int
    .help = "IQR multiplier for determining the threshold value"

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


class RadialProfileSpotFinderThresholdExt:
    """
    Extension to calculate a radial profile threshold. This method calculates
    background value and iqr in 2θ shells, then sets a threshold at a level
    n_iqr above the radial background. As such, it is important to have the
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

        phil = parse(phil_str)
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
        r"""
        Compute the threshold.

        :param image: The image to process
        :param mask: The pixel mask on the image
        :\*\*kwargs: Arbitrary keyword arguments
        :returns: A boolean mask showing foreground/background pixels
        """

        if self.kernel:
            image = convolve(image, self.kernel)

        panel = imageset.get_detector()[i_panel]
        beam = imageset.get_beam()

        # Get 2θ array for the panel or ROI
        two_theta_array = panel.get_two_theta_array(beam.get_s0())
        if region_of_interest:
            x0, x1, y0, y1 = region_of_interest
            two_theta_array = two_theta_array[y0:y1, x0:x1]

        # Convert to 2θ bin selections
        lookup = two_theta_array - flex.min(two_theta_array)
        n_bins = self.params.spotfinder.threshold.radial_profile.n_bins
        multiplier = n_bins / flex.max(lookup + 1e-10)
        lookup *= multiplier  # values now in range [0,n_bins+1)
        lookup = (
            flex.floor(lookup).iround().as_size_t()
        )  # values now in range [0,n_bins-1]

        # Calculate median intensity and IQR within each bin of masked values
        masked_lookup = lookup.select(mask.as_1d())
        masked_image = image.select(mask.as_1d())
        binned_statistics = BinnedStatistics(masked_image, masked_lookup, n_bins)
        med_I = binned_statistics.get_medians()
        iqr = binned_statistics.get_iqrs()

        # Determine the threshold value for each bin. This should be at least
        # 1 quantum greater value than the median to avoid selecting everything
        # in low background cases
        n_iqr = self.params.spotfinder.threshold.radial_profile.n_iqr
        add_level = n_iqr * iqr
        adu = 1 / panel.get_gain()
        add_level.set_selected(add_level <= adu, 2.0 * adu)
        threshold = med_I + add_level

        # Now construct a threshold image
        thresh_im = threshold.select(lookup.as_1d())

        # Peaks are unmasked pixels greater than the threshold
        peaks = image > thresh_im
        return peaks & mask
