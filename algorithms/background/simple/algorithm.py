from __future__ import annotations


class BackgroundAlgorithm:
    """Class to do background subtraction."""

    def __init__(self, experiments, outlier="nsigma", model="constant3d", **kwargs):
        """
        Initialise the algorithm.

        :param experiments: The list of experiments
        :param outlier: The outlier rejection algorithm
        :param model: The background model algorithm
        """
        from dials.algorithms.background.simple import (
            Constant2dModeller,
            Constant3dModeller,
            Creator,
            Linear2dModeller,
            Linear3dModeller,
            MosflmOutlierRejector,
            NormalOutlierRejector,
            NSigmaOutlierRejector,
            TruncatedOutlierRejector,
            TukeyOutlierRejector,
        )

        def select_modeller():
            if model == "constant2d":
                return Constant2dModeller()
            elif model == "constant3d":
                return Constant3dModeller()
            elif model == "linear2d":
                return Linear2dModeller()
            elif model == "linear3d":
                return Linear3dModeller()
            raise RuntimeError(f"Unexpected background model: {model}")

        def select_rejector():
            if outlier == "null":
                return None
            elif outlier == "truncated":
                return TruncatedOutlierRejector(
                    kwargs.get("lower", 0.01), kwargs.get("upper", 0.01)
                )
            elif outlier == "nsigma":
                return NSigmaOutlierRejector(
                    kwargs.get("lower", 3), kwargs.get("upper", 3)
                )
            elif outlier == "normal":
                return NormalOutlierRejector(kwargs.get("min_pixels", 10))
            elif outlier == "plane":
                return MosflmOutlierRejector(
                    kwargs.get("fraction", 1.0), kwargs.get("n_sigma", 4.0)
                )
            elif outlier == "tukey":
                return TukeyOutlierRejector(
                    kwargs.get("lower", 1.5), kwargs.get("upper", 1.5)
                )
            raise RuntimeError(f"Unexpected outlier rejector: {outlier}")

        # Get the minimum number of pixels
        min_pixels = kwargs.get("min_pixels", 10)

        modeller = select_modeller()
        rejector = select_rejector()
        self._creator = Creator(modeller, rejector, min_pixels=min_pixels)

    def compute_background(self, reflections, image_volume=None):
        """
        Compute the background.

        :param reflections: The list of reflections
        """
        from dials.array_family import flex

        # Do the background subtraction
        if image_volume is None:
            reflections["background.mse"] = flex.double(len(reflections))
            reflections["background.dispersion"] = flex.double(len(reflections))
            success = self._creator(
                reflections["shoebox"],
                reflections["background.mse"],
                reflections["background.dispersion"],
            )
            reflections["background.mean"] = reflections["shoebox"].mean_background()
        else:
            success = self._creator(reflections, image_volume)
        reflections.set_flags(~success, reflections.flags.dont_integrate)
        return success


class SimpleBackgroundCalculatorFactory:
    """Class to do background subtraction."""

    @staticmethod
    def create(experiments, outlier="nsigma", model="constant3d", **kwargs):
        """
        Initialise the algorithm.

        :param experiments: The list of experiments
        :param outlier: The outlier rejection algorithm
        :param model: The background model algorithm
        """
        from dials.algorithms.background.simple import (
            Constant2dModeller,
            Constant3dModeller,
            Linear2dModeller,
            Linear3dModeller,
            MosflmOutlierRejector,
            NormalOutlierRejector,
            NSigmaOutlierRejector,
            TruncatedOutlierRejector,
            TukeyOutlierRejector,
        )
        from dials.algorithms.integration.parallel_integrator import (
            SimpleBackgroundCalculator,
        )

        def select_modeller():
            if model == "constant2d":
                return Constant2dModeller()
            elif model == "constant3d":
                return Constant3dModeller()
            elif model == "linear2d":
                return Linear2dModeller()
            elif model == "linear3d":
                return Linear3dModeller()
            raise RuntimeError(f"Unexpected background model: {model}")

        def select_rejector():
            if outlier == "null":
                return None
            elif outlier == "truncated":
                return TruncatedOutlierRejector(
                    kwargs.get("lower", 0.01), kwargs.get("upper", 0.01)
                )
            elif outlier == "nsigma":
                return NSigmaOutlierRejector(
                    kwargs.get("lower", 3), kwargs.get("upper", 3)
                )
            elif outlier == "normal":
                return NormalOutlierRejector(kwargs.get("min_pixels", 10))
            elif outlier == "plane":
                return MosflmOutlierRejector(
                    kwargs.get("fraction", 1.0), kwargs.get("n_sigma", 4.0)
                )
            elif outlier == "tukey":
                return TukeyOutlierRejector(
                    kwargs.get("lower", 1.5), kwargs.get("upper", 1.5)
                )
            raise RuntimeError(f"Unexpected outlier rejector: {outlier}")

        # Get the minimum number of pixels
        min_pixels = kwargs.get("min_pixels", 10)

        modeller = select_modeller()
        rejector = select_rejector()
        return SimpleBackgroundCalculator(modeller, rejector, min_pixels)
