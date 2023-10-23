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
        import os

        from dials.array_family import flex

        # Do the background subtraction
        def print_mask_codes():
            from dials.algorithms.shoebox import MaskCode

            shoebox = reflections["shoebox"]
            nref = 20
            print(f"Shoebox information for the first {nref} reflections.")
            print(
                "sbox size Background BackgroundUsed Foreground Overlapped Strong Valid"
            )
            for i, s in enumerate(shoebox):
                sz = len(s.data)
                b = s.count_mask_values(MaskCode.Background)
                bu = s.count_mask_values(MaskCode.BackgroundUsed)
                f = s.count_mask_values(MaskCode.Foreground)
                o = s.count_mask_values(MaskCode.Overlapped)
                st = s.count_mask_values(MaskCode.Strong)
                v = s.count_mask_values(MaskCode.Valid)
                print(
                    f"{i:4d} {sz:4d}       {b:4d}           {bu:4d}       {f:4d}       {o:4d}   {st:4d}  {v:4d}"
                )
                if i == nref:
                    break

        if "COMPUTE_BACKGROUND" in os.environ:
            print(
                "In BackgroundAlgorithm.compute_background *before* background calculation"
            )
            print_mask_codes()

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

        if "COMPUTE_BACKGROUND" in os.environ:
            print(
                f"Background calculation success for {success.count(True)} of {len(reflections)} reflections"
            )
            print(
                "In BackgroundAlgorithm.compute_background *after* background calculation"
            )
            print_mask_codes()

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
