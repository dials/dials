from __future__ import absolute_import, division, print_function


class ComputeDetectorNoise(object):
    """Calculate the detector noise. Estimate this by calculating the mean
    of the corner pixels in a series of images and return a constant."""

    def __init__(self):
        """Initialise the algorithm."""
        from scitbx.array_family import flex

        self._pixels = flex.int()

    def add(self, image):
        """Add another image to the calculation.

        Params:
            image The image to use
        """
        # Get size of the image
        height, width = image.all()

        # Add pixels to the list
        self._pixels.append(image[0, 0])
        self._pixels.append(image[0, width - 1])
        self._pixels.append(image[height - 1, 0])
        self._pixels.append(image[height - 1, width - 1])

    def pixels(self):
        """Get the pixels used in the calculation.

        Returns:
            The array of pixel intensities
        """
        return self._pixels

    def compute(self):
        """Compute the noise.

        Returns:
            The calculated detector noise
        """
        from scitbx.array_family import flex

        return flex.mean(self._pixels.as_double())
