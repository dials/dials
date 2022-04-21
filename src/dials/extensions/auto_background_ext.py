from __future__ import annotations

from dials.extensions.glm_background_ext import GLMBackgroundExt
from dials.extensions.simple_background_ext import SimpleBackgroundExt


class AutoBackgroundExt:
    """An extension class to automatically choose a background algorithm.
    If experiments all have counting detectors then the glm algorithm is
    chosen. Otherwise, the simple algorithm is chosen."""

    name = "Auto"

    default = True

    def __init__(self, params, experiments):
        """
        Initialise the algorithm.

        :param params: The input parameters
        :param experiments: The list of experiments
        """

        # Guess at whether we have a counting detector, based on gains and
        # detector types (https://github.com/dials/dials/issues/706)
        detectors = experiments.detectors()
        gains = [p.get_gain() for d in detectors for p in d]
        pnl_types = [p.get_type() for d in detectors for p in d]
        counting_detector = True
        if [g == 1.0 for g in gains].count(False) != 0:
            counting_detector = False
        if ["PAD" in t for t in pnl_types].count(False) != 0:
            counting_detector = False

        if counting_detector:
            ext = GLMBackgroundExt(params, experiments)
        else:
            ext = SimpleBackgroundExt(params, experiments)

        self._algorithm = ext._algorithm

    def compute_background(self, reflections, image_volume=None):
        """
        Compute the background.

        :param reflections: The list of reflections
        """
        return self._algorithm.compute_background(
            reflections, image_volume=image_volume
        )
