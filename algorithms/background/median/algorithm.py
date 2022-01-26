from __future__ import annotations


class BackgroundAlgorithm:
    """Class to do background subtraction."""

    def __init__(self, experiments):
        """
        Initialise the algorithm.

        :param experiments: The list of experiments
        """
        pass

    def compute_background(self, reflections, image_volume=None):
        """
        Compute the background.

        :param reflections: The list of reflections
        """
        from dials.algorithms.background.median import create
        from dials.array_family import flex

        # Do the background subtraction
        if image_volume is None:
            success = create(reflections["shoebox"])
            reflections["background.mean"] = flex.double(
                [sbox.background[0] for sbox in reflections["shoebox"]]
            )
        else:
            success = create(reflections, image_volume)
        reflections.set_flags(~success, reflections.flags.dont_integrate)
        return success
