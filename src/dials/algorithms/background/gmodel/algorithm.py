from __future__ import annotations

import pickle


class ModelCache:
    """
    A class to cache the model
    """

    def __init__(self):
        """
        Create a model dictionary
        """
        self.model = {}

    def get(self, name):
        """
        Get the model
        """
        if name is None:
            raise RuntimeError("Model is not specified")
        try:
            model = self.model[name]
        except KeyError:
            with open(name, "rb") as infile:
                model = pickle.load(infile)
                self.model[name] = model
        return model


# Instance of the model cache
global_model_cache = ModelCache()


class BackgroundAlgorithm:
    """Class to do background subtraction."""

    def __init__(
        self,
        experiments,
        model=None,
        robust=False,
        tuning_constant=1.345,
        min_pixels=10,
    ):
        """
        Initialise the algorithm.

        :param experiments: The list of experiments
        :param model: The background model
        :param robust: Use the robust background algorithm
        :param tuning_constant: The robust tuning constant
        """
        from dials.algorithms.background.gmodel import Creator

        # Get the model
        model = global_model_cache.get(model)

        # Create the background creator
        self._create = Creator(model=model, robust=robust, min_pixels=min_pixels)

    def compute_background(self, reflections, image_volume=None):
        """
        Compute the background.

        :param reflections: The list of reflections
        """
        # Do the background subtraction
        if image_volume is None:
            success = self._create(reflections)
            reflections["background.mean"] = reflections[
                "shoebox"
            ].mean_modelled_background()
        else:
            success = self._create(reflections, image_volume)
        reflections.set_flags(~success, reflections.flags.dont_integrate)
        return success


class GModelBackgroundCalculatorFactory:
    """Class to do background subtraction."""

    @staticmethod
    def create(experiments, model=None, robust=False, min_pixels=10):
        """
        Initialise the algorithm.

        :param experiments: The list of experiments
        :param model: The background model
        :param robust: Use the robust background algorithm
        :param tuning_constant: The robust tuning constant
        """
        from dials.algorithms.integration.parallel_integrator import (
            GModelBackgroundCalculator,
        )

        # Get the model
        model = global_model_cache.get(model)

        # Create the background creator
        return GModelBackgroundCalculator(
            model=model, robust=robust, min_pixels=min_pixels
        )
