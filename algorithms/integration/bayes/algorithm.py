from __future__ import annotations


class IntegrationAlgorithm:
    """A class to perform bayesian integration"""

    def __init__(self, **kwargs):
        pass

    def __call__(self, reflections, image_volume=None):
        """Process the reflections.

        :param reflections: The reflections to integrate
        :return: The list of integrated reflections
        """
        # Integrate and return the reflections
        if image_volume is None:
            intensity = reflections["shoebox"].bayesian_summation_intensity()
        else:
            raise RuntimeError("Image volume not supported at the moment")
        reflections["intensity.sum.value"] = intensity.observed_value()
        reflections["intensity.sum.variance"] = intensity.observed_variance()
        reflections["background.sum.value"] = intensity.background_value()
        reflections["background.sum.variance"] = intensity.background_variance()
        success = intensity.observed_success()
        reflections.set_flags(success, reflections.flags.integrated_sum)
        return success
