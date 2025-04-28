"""Define a class that searches for beam position using midpoint method"""

from __future__ import annotations

import numpy as np

from dials.algorithms.beam_position.helper_functions import normalize
from dials.algorithms.beam_position.project_profile import project


class InversionMethodSolver:
    def __init__(self, image, params, axis="x"):
        background_cutoff = params.projection.inversion.background_cutoff
        threshold = params.projection.inversion.bad_pixel_threshold
        self.axis = axis

        if axis == "x":
            exclude_range = params.projection.exclude_pixel_range_y
        elif axis == "y":
            exclude_range = params.projection.exclude_pixel_range_x
        else:
            raise ValueError(f"Unknown axis: {axis}")

        clean_image = np.array(image)

        if threshold:
            clean_image[clean_image > threshold] = 0.0
        if background_cutoff:
            clean_image[clean_image < background_cutoff] = 0.0

        profile_max, pmax, pmin = project(
            clean_image,
            axis=axis,
            method="max",
            exclude_range=exclude_range,
            convolution_width=1,
            normalize=False,
        )

        self.params = params
        self.profile_max = profile_max
        self.max_value = pmax
        self.min_value = pmin

    def find_beam_position(self):
        n = len(self.profile_max)
        guess_position = self.params.projection.inversion.guess_position
        window_width = self.params.projection.inversion.inversion_window_width

        if guess_position:
            if len(guess_position) != 2:
                msg = "guess_position must be a tuple of two ints"
                raise ValueError(msg)
            guess_x, guess_y = guess_position
            if self.axis == "x":
                center = guess_position[0]
            else:
                center = guess_position[1]
        else:
            center = int(n / 2)

        min_ind = center - window_width
        if min_ind < 0:
            min_ind = 0
        max_ind = center + window_width

        if max_ind > n:
            max_ind = n - 1

        indices = np.arange(min_ind, max_ind, 1)

        correlations = 0 * self.profile_max

        for i in indices:
            correlations[i] = invert_and_correlate(self.profile_max, i)

        scale = 0.8 * self.profile_max.max()  # Scale for plotting
        correlations = scale * normalize(correlations)

        self.beam_position = correlations.argmax()
        self.correlations = correlations

        return float(self.beam_position)

    def plot(self, figure):
        n = len(self.profile_max)
        indices = np.arange(n)

        if self.axis == "x":
            ax = figure.axis_x
            ax.plot(indices, self.profile_max, lw=1, c="gray", label="max. projection")

            ax.plot(indices, self.correlations, lw=1, c="C2", label="invert and sum")

            ax.text(
                0.01,
                0.95,
                "method: inversion",
                va="top",
                ha="left",
                transform=ax.transAxes,
                fontsize=7,
            )
            label = f"I_min, I_max = ({self.min_value:.1f}, {self.max_value:.1f})"
            ax.text(
                0.01,
                0.78,
                label,
                va="top",
                ha="left",
                c="gray",
                transform=ax.transAxes,
                fontsize=5,
            )
            ax.legend(
                loc=(0.6, 0.7),
                labelspacing=0.5,
                borderpad=0,
                columnspacing=3.5,
                handletextpad=0.4,
                fontsize=7,
                handlelength=2.0,
                handleheight=0.7,
                frameon=False,
            )
            ax.axvline(self.beam_position, c="C3", lw=1)
            ax.tick_params(axis="y", colors="gray")

        elif self.axis == "y":
            ax = figure.axis_y
            ax.plot(self.profile_max, indices, lw=1, c="gray", label="max. projection")

            ax.plot(self.correlations, indices, lw=1, c="C2", label="invert and sum")

            ax.text(
                0.95,
                0.99,
                "method: inversion",
                va="top",
                ha="right",
                transform=ax.transAxes,
                rotation=-90,
                fontsize=7,
            )
            label = f"I_min, I_max = ({self.min_value:.1f}, {self.max_value:.1f})"
            ax.text(
                0.78,
                0.99,
                label,
                va="top",
                ha="right",
                transform=ax.transAxes,
                rotation=-90,
                fontsize=5,
            )
            ax.legend(
                loc=(0.02, 0.1),
                labelspacing=0.5,
                borderpad=0,
                columnspacing=3.5,
                handletextpad=0.4,
                fontsize=7,
                handlelength=1.5,
                handleheight=0.7,
                frameon=False,
            )
            ax.axhline(self.beam_position, c="C3", lw=1)
            ax.tick_params(axis="x", colors="gray")
        else:
            raise ValueError(f"Unknown axis: {self.axis}")

        # Set matplotlib legend
        # handles, labels = ax.get_legend_handles_labels()


def invert_and_correlate(x, index):
    """
    Given an 1D array x, compute inverted array inv_x (inversion around an
    element with an index 'index') and return a sum of x * inv_x.
    """

    inv_x = np.zeros(len(x))
    n = len(x)
    half = int(n / 2)

    # Compute the inverted 1D array
    if index <= half:
        left = x[0:index]
        right = x[index : 2 * index]
        inv_x[0:index] = right[::-1]
        inv_x[index : 2 * index] = left[::-1]
    else:
        right = x[index:]
        width = len(right)
        left = x[index - width : index]
        inv_x[index - width : index] = right[::-1]
        inv_x[index:] = left[::-1]

    return np.mean(x * inv_x)  # Could try with mean instead of sum
