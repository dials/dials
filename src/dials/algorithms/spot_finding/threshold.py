from __future__ import annotations


class ThresholdStrategy:
    """
    Base class for spot finder threshold strategies.
    """

    def __init__(self, **kwargs):
        """
        Initialise with key word arguments.
        """
        pass

    def __call__(self, image):
        """
        Threshold the image.
        """
        raise RuntimeError("Overload Me!")


class DispersionThresholdStrategy(ThresholdStrategy):
    """
    A class implementing a 'gain' threshold.
    """

    def __init__(self, **kwargs):
        """
        Set the threshold algorithm up
        """

        # Initialise the base class
        ThresholdStrategy.__init__(self, **kwargs)

        # Get the parameters
        self._kernel_size = kwargs.get("kernel_size", (3, 3))
        self._gain = kwargs.get("gain")
        self._n_sigma_b = kwargs.get("n_sigma_b", 6)
        self._n_sigma_s = kwargs.get("n_sigma_s", 3)
        self._min_count = kwargs.get("min_count", 2)
        self._threshold = kwargs.get("global_threshold", 0)

        # Save the constant gain
        self._gain_map = None

        # Create a buffer
        self.algorithm = {}

    def __call__(self, image, mask):
        """
        Call the thresholding function

        :param image: The image to process
        :param mask: The mask to use
        :return: The thresholded image
        """
        from dials.algorithms.image import threshold
        from dials.array_family import flex

        # Initialise the algorithm
        try:
            algorithm = self.algorithm[image.all()]
        except Exception:
            algorithm = threshold.DispersionThreshold(
                image.all(),
                self._kernel_size,
                self._n_sigma_b,
                self._n_sigma_s,
                self._threshold,
                self._min_count,
            )
            self.algorithm[image.all()] = algorithm

        # Set the gain
        if self._gain is not None:
            assert self._gain > 0
            self._gain_map = flex.double(image.accessor(), self._gain)
            self._gain = None

        # Compute the threshold
        result = flex.bool(flex.grid(image.all()))
        if self._gain_map:
            algorithm(image, mask, self._gain_map, result)
        else:
            algorithm(image, mask, result)

        # Return the result
        return result


class DispersionExtendedThresholdStrategy(ThresholdStrategy):
    """
    A class implementing a 'gain' threshold.
    """

    def __init__(self, **kwargs):
        """
        Set the threshold algorithm up
        """

        # Initialise the base class
        ThresholdStrategy.__init__(self, **kwargs)

        # Get the parameters
        self._kernel_size = kwargs.get("kernel_size", (3, 3))
        self._gain = kwargs.get("gain")
        self._n_sigma_b = kwargs.get("n_sigma_b", 6)
        self._n_sigma_s = kwargs.get("n_sigma_s", 3)
        self._min_count = kwargs.get("min_count", 2)
        self._threshold = kwargs.get("global_threshold", 0)

        # Save the constant gain
        self._gain_map = None

        # Create a buffer
        self.algorithm = {}

    def __call__(self, image, mask):
        """
        Call the thresholding function

        :param image: The image to process
        :param mask: The mask to use
        :return: The thresholded image
        """
        from dials.algorithms.image import threshold
        from dials.array_family import flex

        # Initialise the algorithm
        try:
            algorithm = self.algorithm[image.all()]
        except Exception:
            algorithm = threshold.DispersionExtendedThreshold(
                image.all(),
                self._kernel_size,
                self._n_sigma_b,
                self._n_sigma_s,
                self._threshold,
                self._min_count,
            )
            self.algorithm[image.all()] = algorithm

        # Set the gain
        if self._gain is not None:
            assert self._gain > 0
            self._gain_map = flex.double(image.accessor(), self._gain)
            self._gain = None

        # Compute the threshold
        result = flex.bool(flex.grid(image.all()))
        if self._gain_map:
            algorithm(image, mask, self._gain_map, result)
        else:
            algorithm(image, mask, result)

        # Return the result
        return result
