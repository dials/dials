from __future__ import absolute_import, division, print_function


class IntegrationAlgorithm(object):
    """A class to perform 3D summation integration"""

    def __init__(self, **kwargs):
        pass

    def __call__(self, reflections, image_volume=None):
        """Process the reflections.

        :param reflections: The reflections to integrate
        :return: The list of integrated reflections
        """
        from dials.algorithms.integration.sum import sum_image_volume

        # Integrate and return the reflections
        if image_volume is None:
            intensity = reflections["shoebox"].summed_intensity()
        else:
            intensity = sum_image_volume(reflections, image_volume)
        reflections["intensity.sum.value"] = intensity.observed_value()
        reflections["intensity.sum.variance"] = intensity.observed_variance()
        reflections["background.sum.value"] = intensity.background_value()
        reflections["background.sum.variance"] = intensity.background_variance()
        success = intensity.observed_success()
        reflections.set_flags(success, reflections.flags.integrated_sum)
        return success
