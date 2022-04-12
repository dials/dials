from __future__ import annotations


class SimpleCentroidExt:
    """An extension class implementing a simple centroid algorithm."""

    name = "simple"

    default = True

    def __init__(self, params, experiments):
        """Initialise the algorithm.

        :param params: The input phil parameters
        :param experiments: The experiment list
        """
        self.experiments = experiments

    def compute_centroid(self, reflections, image_volume=None):
        """
        Compute the centroid.

        :param reflections: The list of reflections
        """
        import dials.algorithms.centroid.simple

        return dials.algorithms.centroid.simple.centroid(
            self.experiments, reflections, image_volume=image_volume
        )
