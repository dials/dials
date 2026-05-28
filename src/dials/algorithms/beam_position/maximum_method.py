"""Define a class that searches for beam position using midpoint method"""

from __future__ import annotations

from collections import namedtuple

import numpy as np
from matplotlib.ticker import MultipleLocator

from dials.algorithms.beam_position.project_profile import project


class MaximumMethodSolver:
    def __init__(self, image, params, axis="x"):
        cw = params.projection.maximum.convolution_width
        threshold = params.projection.maximum.bad_pixel_threshold
        nc = params.projection.maximum.n_convolutions
        self.axis = axis

        if axis == "x":
            exclude_range = params.projection.exclude_pixel_range_y
        elif axis == "y":
            exclude_range = params.projection.exclude_pixel_range_x
        else:
            raise ValueError(f"Unknown axis: {axis}")

        # Clean image
        clean_image = np.array(image)
        min_value = clean_image.min()

        if threshold:
            clean_image[clean_image > threshold] = min_value

        profile_max, pmax, pmin = project(
            clean_image,
            axis=axis,
            method="max",
            exclude_range=exclude_range,
            convolution_width=1,
        )

        profile_mean, pmax1, pmin1 = project(
            clean_image,
            axis=axis,
            method="average",
            exclude_range=exclude_range,
            convolution_width=cw,
            n_convolutions=nc,
        )

        self.params = params
        self.profile_max = profile_max
        self.profile_mean = profile_mean
        self.max_value = pmax
        self.min_value = pmin
        self.max_from_mean = pmax1
        self.min_from_mean = pmin1
        self.bin_width = params.projection.maximum.bin_width
        self.bin_step = params.projection.maximum.bin_step
        if self.bin_width < self.bin_step:
            msg = (
                f"maximum.bin_width = {self.bin_width} is smaller then the "
                f"maximum.bin_step = {self.bin_step}."
            )
            raise ValueError(msg)
        self.n_convolutions = params.projection.maximum.n_convolutions

    def find_beam_position(self):
        Bin = namedtuple("Bin", ["value", "start", "end"])

        n = len(self.profile_max)

        if self.bin_width > n:
            msg = "Error in MaximumMethodSolver: bin_width="
            msg += f"{self.bin_width} larger than image dimension {n}"
            raise ValueError(msg)

        n_end = ((n - self.bin_width) // self.bin_step) * self.bin_step

        bins = []

        for start in np.arange(0, n_end, self.bin_step):
            value = self.profile_mean[start : start + self.bin_width].sum()
            bins.append(Bin(value, start, start + self.bin_width))

        max_bin = sorted(bins, key=lambda x: x.value, reverse=True)[0]

        selected = np.array(self.profile_max)
        min_value = selected.min()
        selected[0 : max_bin.start] = min_value
        selected[max_bin.end :] = min_value

        beam_position = np.argmax(selected)

        self.beam_position = beam_position
        self.bin_start = max_bin.start
        self.bin_end = max_bin.end

        return beam_position

    def plot(self, figure):
        indices = np.arange(0, len(self.profile_max))

        if self.axis == "x":
            ax = figure.axis_x
            ax.axvspan(self.bin_start, self.bin_end, color="#BEBEBE", alpha=0.5, lw=0)

            ax.text(
                0.01,
                0.78,
                f"I_min, I_max = ({self.min_value:.1f}, {self.max_value:.1f})",
                ha="left",
                transform=ax.transAxes,
                c="gray",
                fontsize=5,
            )
            ax.text(
                0.01,
                0.7,
                f"I_min, I_max = ({self.min_from_mean:.1f}, {self.max_from_mean:.1f})",
                ha="left",
                transform=ax.transAxes,
                c="C2",
                fontsize=5,
            )

            ax.plot(indices, self.profile_max, lw=1, c="gray", label="max. projection")
            ax.plot(indices, self.profile_mean, lw=1, c="C2", label="avg. projection")

            ax.text(
                0.01,
                0.95,
                "method: maximum",
                va="top",
                ha="left",
                transform=ax.transAxes,
                fontsize=7,
            )

            ax.legend(
                loc=(0.6, 0.7),
                labelspacing=0.5,
                borderpad=0,
                columnspacing=3.5,
                handletextpad=0.4,
                fontsize=7,
                handlelength=1.5,
                handleheight=0.7,
                frameon=False,
            )

            ax.axvline(self.beam_position, c="C3", lw=1)

            ax.set_ylim(-0.07, 1.2)
            ax.set_yticks([0, 0.5, 1])

            mloc = MultipleLocator(0.1)
            ax.yaxis.set_minor_locator(mloc)
            ax.tick_params(axis="y", colors="black")

        elif self.axis == "y":
            ax = figure.axis_y
            ax.axhspan(self.bin_start, self.bin_end, color="#BEBEBE", alpha=0.5, lw=0)
            ax.plot(self.profile_max, indices, lw=1, c="gray", label="max. projection")
            ax.plot(self.profile_mean, indices, lw=1, c="C2", label="avg. projection")
            ax.text(
                0.95,
                0.99,
                "method: maximum",
                va="top",
                ha="right",
                transform=ax.transAxes,
                rotation=-90,
                fontsize=7,
            )

            ax.text(
                0.0,
                1.06,
                f"I_min, I_max = ({self.min_value:.1f}, {self.max_value:.1f})",
                va="top",
                transform=ax.transAxes,
                c="gray",
                fontsize=5,
            )
            ax.text(
                0.0,
                1.03,
                f"I_min, I_max = ({self.min_from_mean:.1f}, {self.max_from_mean:.1f})",
                va="top",
                transform=ax.transAxes,
                c="C2",
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

            ax.set_xlim(-0.07, 1.2)
            ax.set_xticks([0, 0.5, 1])

            mloc = MultipleLocator(0.1)
            ax.xaxis.set_minor_locator(mloc)

            ax.tick_params(axis="x", colors="black")
        else:
            raise ValueError(f"Unknown axis: {self.axis}")
