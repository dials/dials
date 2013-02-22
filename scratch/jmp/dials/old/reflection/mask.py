
class ReflectionMask:
    """A class to represent the reflection mask."""

    def __init__(self, sigma_divergence, sigma_mosaicity, n_sigma=10):
        """Initialise the mask parameters.

        Args:
            sigma_divergence: The standard deviation of the beam divergence
            sigma_mosaicity: The standard deviation of the mosaicity
            n_sigma: The size of the mask

        """
        self.n_sigma = n_sigma
        self.sigma_divergence = sigma_divergence
        self.sigma_mosaicity = sigma_mosaicity
        self.delta_divergence = self.n_sigma * self.sigma_divergence
        self.delta_mosaicity = self.n_sigma * self.sigma_mosaicity

    def is_in_mask(self, ep1, ep2, ep3):
        """Check if the point is in the mask or not.

        TODO:
            For the dataset I've been using, excluding along e3 sometimes means
            1 or 0 data frames are used!

        Args:
            ep1: The e1 coordinate
            ep2: The e2 coordinate
            ep3: The e3 coordinate

        Returns:
            True/False

        """
        return (abs(ep1) <= 0.5 * self.delta_divergence and
                abs(ep2) <= 0.5 * self.delta_divergence)# and
                #abs(ep3) <= 0.5 * self.delta_mosaicity)
