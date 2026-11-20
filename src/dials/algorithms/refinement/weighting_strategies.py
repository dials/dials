"""Contains classes used to provide weighting schemes as strategies for
ReflectionManagers."""

from __future__ import annotations

import numpy as np

from dxtbx.model.tof_helpers import frame_to_tof_interpolator

from dials.algorithms.refinement import DialsRefineConfigError
from dials.array_family import flex


class StatisticalWeightingStrategy:
    """Defines a single method that provides a ReflectionManager with a strategy
    for calculating weights for refinement"""

    @staticmethod
    def calculate_weights(reflections):
        """set 'statistical weights', that is w(x) = 1/var(x)"""

        weights = (reflections["xyzobs.mm.variance"]).deep_copy()
        parts = weights.parts()
        for w in parts:
            sel = w > 0.0
            w.set_selected(sel, 1.0 / w.select(sel))
        reflections["xyzobs.mm.weights"] = flex.vec3_double(*parts)
        indexed = reflections.select(reflections.get_flags(reflections.flags.indexed))
        if any(indexed["xyzobs.mm.weights"].norms() == 0.0):
            raise DialsRefineConfigError(
                "Cannot set statistical weights as some indexed reflections have observed variances equal to zero"
            )

        return reflections


class StillsWeightingStrategy(StatisticalWeightingStrategy):
    """Defines a single method that provides a ReflectionManager with a strategy
    for calculating weights for refinement. This version uses statistical weights
    for X and Y and a fixed constant for the delta Psi part, defaulting to 1000000"""

    def __init__(self, delpsi_constant=1000000):
        self._delpsi_constant = delpsi_constant

    def calculate_weights(self, reflections):
        """Include weights for DeltaPsi"""

        # call parent class method to set X and Y weights
        reflections = super().calculate_weights(reflections)

        reflections["delpsical.weights"] = flex.double(
            len(reflections), self._delpsi_constant
        )

        return reflections


class ExternalDelPsiWeightingStrategy(StatisticalWeightingStrategy):
    """Defines a single method that provides a ReflectionManager with a strategy
    for calculating weights for stills refinement. This version uses statistical
    weights for X and Y and assume that the Delta Psi part is already provided in
    the reflection table"""

    def calculate_weights(self, reflections):
        """Statistical weights for X, Y. Weights for DeltaPsi must be already
        provided in the reflection table"""

        # call parent class method to set X and Y weights
        reflections = super().calculate_weights(reflections)

        if "delpsical.weights" not in reflections:
            raise DialsRefineConfigError(
                'The key "delpsical.weights" is expected within the input reflections'
            )

        return reflections


class ConstantWeightingStrategy:
    def __init__(self, wx, wy, wz, stills=False):
        self._wx = wx
        self._wy = wy
        self._wz = wz
        self._stills = stills

    def calculate_weights(self, reflections):
        """Set weights to constant terms. If stills, the z weights are
        the 'delpsical.weights' attribute of the reflection table. Otherwise, use
        the usual 'xyzobs.mm.weights'"""

        wx = flex.double(len(reflections), self._wx)
        wy = flex.double(len(reflections), self._wy)
        wz = flex.double(len(reflections), self._wz)
        if self._stills:
            null = flex.double(len(reflections), 0)
            reflections["xyzobs.mm.weights"] = flex.vec3_double(wx, wy, null)
            reflections["delpsical.weights"] = wz
        else:
            reflections["xyzobs.mm.weights"] = flex.vec3_double(wx, wy, wz)

        return reflections


class LaueStatisticalWeightingStrategy(StatisticalWeightingStrategy):
    """
    Variance in z estimated from sqrt(x^2+y^2)
    """

    def __init__(
        self,
        wavelength_weight: float = 1e7,
    ):
        self._wavelength_weight = wavelength_weight

    def calculate_weights(self, reflections):
        reflections = super().calculate_weights(reflections)

        wx, wy, _ = reflections["xyzobs.mm.weights"].parts()
        wz = flex.sqrt(wx * wx + wy * wy) * self._wavelength_weight
        reflections["xyzobs.mm.weights"] = flex.vec3_double(wx, wy, wz)

        return reflections


class LaueMixedWeightingStrategy(StatisticalWeightingStrategy):
    """
    Use statistical weighting for x and y, and constant weighting for z
    """

    def __init__(
        self,
        wavelength_weight: float = 1e7,
    ):
        self._wavelength_weight = wavelength_weight

    def calculate_weights(self, reflections):
        reflections = super().calculate_weights(reflections)

        wx, wy, wz = reflections["xyzobs.mm.weights"].parts()
        wz = wz * 0 + self._wavelength_weight
        reflections["xyzobs.mm.weights"] = flex.vec3_double(wx, wy, wz)

        return reflections


class TOFStatisticalWeightingStrategy:
    """Weights for time-of-flight Laue refinement (X, Y (mm) and ToF (usec)).

    X and Y use statistical weights, w = 1/var, from the centroid variances.

    The observed ToF is the centre of the maximum ToF bin of the spot, and
    the ToF uncertainty is taken as one bin width
    at the observed position, where sigma_t = dt/dframe, and w = 1/sigma_t^2."""

    def __init__(self, experiments):
        self._experiments = experiments

    def _tof_bin_widths(self, reflections):
        """ToF bin width (usec) at the observed frame of each reflection"""

        bin_widths = flex.double(len(reflections), 0.0)
        frames_obs = reflections["xyzobs.px.value"].parts()[2]
        for iexp, expt in enumerate(self._experiments):
            sel = reflections["id"] == iexp
            if sel.count(True) == 0:
                continue
            tof = expt.scan.get_property("time_of_flight")  # (usec)
            frames = list(range(len(tof)))
            frame_to_tof = frame_to_tof_interpolator(frames, tof)

            # central difference of +/- half a frame, clipped to the frame range
            z = frames_obs.select(sel).as_numpy_array()
            lo = np.clip(z - 0.5, frames[0], frames[-1])
            hi = np.clip(z + 0.5, frames[0], frames[-1])
            dtof_dframe = (frame_to_tof(hi) - frame_to_tof(lo)) / (hi - lo)
            bin_widths.set_selected(sel, flex.double(dtof_dframe))
        return bin_widths

    def calculate_weights(self, reflections):
        var_x, var_y, _ = reflections["xyzobs.mm.variance"].parts()
        wx = var_x.deep_copy()
        wy = var_y.deep_copy()
        for w in (wx, wy):
            sel = w > 0.0
            w.set_selected(sel, 1.0 / w.select(sel))

        sigma_t = self._tof_bin_widths(reflections)
        wt = flex.double(len(reflections), 0.0)
        sel = sigma_t > 0.0
        wt.set_selected(sel, 1.0 / (sigma_t.select(sel) ** 2))

        reflections["xyzobs.mm.weights"] = flex.vec3_double(wx, wy, wt)

        indexed = reflections.select(reflections.get_flags(reflections.flags.indexed))
        if any(indexed["xyzobs.mm.weights"].norms() == 0.0):
            raise DialsRefineConfigError(
                "Cannot set statistical weights as some indexed reflections have "
                "zero weights in X, Y and ToF"
            )

        return reflections


class TOFConstantWeightingStrategy:
    """Constant weights for time-of-flight (ToF) Laue refinement, for the X, Y
    (mm) and ToF (usec) residuals"""

    def __init__(self, wx, wy, wtof):
        self._wx = wx
        self._wy = wy
        self._wtof = wtof

    def calculate_weights(self, reflections):
        n = len(reflections)
        reflections["xyzobs.mm.weights"] = flex.vec3_double(
            flex.double(n, self._wx),
            flex.double(n, self._wy),
            flex.double(n, self._wtof),
        )
        return reflections
