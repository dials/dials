from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger("dials.extensions.dispersion_spotfinder_threshold_ext")


class DispersionSpotFinderThresholdExt(object):
    """ Extensions to do dispersion threshold. """

    name = "dispersion"

    default = True

    @classmethod
    def phil(cls):
        from libtbx.phil import parse

        phil = parse(
            """
      gain = None
        .type = float(value_min=0.0)
        .help = "Use a flat gain map for the entire detector to act as a"
                "multiplier for the gain set by the format. Cannot be used"
                "in conjunction with lookup.gain_map parameter."

      kernel_size = 3 3
        .help = "The size of the local area around the spot in which"
                "to calculate the mean and variance. The kernel is"
                "given as a box of size (2 * nx + 1, 2 * ny + 1) centred"
                "at the pixel."
        .type = ints(size=2)
        .expert_level = 1

      sigma_background = 6
        .help = "The number of standard deviations of the index of dispersion"
                "(variance / mean) in the local area below"
                "which the pixel will be classified as background."
        .type = float
        .expert_level = 1

      sigma_strong = 3
        .help = "The number of standard deviations above the mean in the"
                "local area above which the pixel will be classified as"
                "strong."
        .type = float
        .expert_level = 1

      min_local = 2
        .help = "The minimum number of pixels under the image processing kernel"
                "that are need to do the thresholding operation. Setting the"
                "value between 2 and the total number of pixels under the"
                "kernel will force the algorithm to use that number as the"
                "minimum. If the value is less than or equal to zero, then"
                "the algorithm will use all pixels under the kernel. In"
                "effect this will add a border of pixels which are always"
                "classed as background around the edge of the image and around"
                "any masked out pixels."
        .type = int
        .expert_level = 1

      global_threshold = 0
        .type = float
        .help = "The global threshold value. Consider all pixels less than this"
                "value to be part of the background."
    """
        )
        return phil

    def __init__(self, params):
        """
        Initialise the algorithm.

        :param params: The input parameters

        """
        self.params = params

    def compute_threshold(self, image, mask):
        """
        Compute the threshold.

        :param image: The image to process
        :param mask: The pixel mask on the image
        :returns: A boolean mask showing foreground/background pixels

        """

        import libtbx

        params = self.params
        if params.spotfinder.threshold.dispersion.global_threshold is libtbx.Auto:
            params.spotfinder.threshold.dispersion.global_threshold = int(
                estimate_global_threshold(image, mask)
            )
            logger.info(
                "Setting global_threshold: %i"
                % (params.spotfinder.threshold.dispersion.global_threshold)
            )

        from dials.algorithms.spot_finding.threshold import DispersionThresholdStrategy

        self._algorithm = DispersionThresholdStrategy(
            kernel_size=params.spotfinder.threshold.dispersion.kernel_size,
            gain=params.spotfinder.threshold.dispersion.gain,
            mask=params.spotfinder.lookup.mask,
            n_sigma_b=params.spotfinder.threshold.dispersion.sigma_background,
            n_sigma_s=params.spotfinder.threshold.dispersion.sigma_strong,
            min_count=params.spotfinder.threshold.dispersion.min_local,
            global_threshold=params.spotfinder.threshold.dispersion.global_threshold,
        )

        return self._algorithm(image, mask)


def estimate_global_threshold(image, mask=None, plot=False):

    from scitbx.array_family import flex
    from scitbx import matrix

    n_above_threshold = flex.size_t()
    threshold = flex.double()
    for i in range(1, 20):
        g = 1.5 ** i
        g = int(g)
        n_above_threshold.append((image > g).count(True))
        threshold.append(g)

    # Find the elbow point of the curve, in the same manner as that used by
    # distl spotfinder for resolution method 1 (Zhang et al 2006).
    # See also dials/algorithms/spot_finding/per_image_analysis.py

    x = threshold.as_double()
    y = n_above_threshold.as_double()
    slopes = (y[-1] - y[:-1]) / (x[-1] - x[:-1])
    p_m = flex.min_index(slopes)

    x1 = matrix.col((x[p_m], y[p_m]))
    x2 = matrix.col((x[-1], y[-1]))

    gaps = flex.double()
    v = matrix.col(((x2[1] - x1[1]), -(x2[0] - x1[0]))).normalize()

    for i in range(p_m, len(x)):
        x0 = matrix.col((x[i], y[i]))
        r = x1 - x0
        g = abs(v.dot(r))
        gaps.append(g)

    p_g = flex.max_index(gaps)

    x_g_ = x[p_g + p_m]
    y_g_ = y[p_g + p_m]

    # more conservative, choose point 2 left of the elbow point
    x_g = x[p_g + p_m - 2]
    y_g = y[p_g + p_m - 2]

    if plot:
        from matplotlib import pyplot

        pyplot.figure(figsize=(16, 12))
        pyplot.scatter(threshold, n_above_threshold, marker="+")
        # for i in range(len(threshold)-1):
        #  pyplot.plot([threshold[i], threshold[-1]],
        #              [n_above_threshold[i], n_above_threshold[-1]])
        # for i in range(1, len(threshold)):
        #  pyplot.plot([threshold[0], threshold[i]],
        #              [n_above_threshold[0], n_above_threshold[i]])
        pyplot.plot([x_g, y_g], pyplot.ylim())
        pyplot.plot(
            [threshold[p_m], threshold[-1]],
            [n_above_threshold[p_m], n_above_threshold[-1]],
        )
        pyplot.plot([x_g_, threshold[-1]], [y_g_, n_above_threshold[-1]])
        pyplot.xlabel("Threshold")
        pyplot.ylabel("Number of pixels above threshold")
        pyplot.savefig("global_threshold.png")
        pyplot.clf()

    return x_g
