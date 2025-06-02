from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from matplotlib.patches import Circle

from dials.algorithms.beam_position.helper_functions import remove_pixels_by_intensity
from dials.algorithms.beam_position.project_profile import convert_range_into_spans


class Figure:
    def __init__(self, filename="fig.png"):
        fig = plt.figure(figsize=(5, 5))
        gs = gridspec.GridSpec(
            2,
            2,
            top=0.90,
            bottom=0.09,
            left=0.15,
            right=0.96,
            wspace=0.04,
            hspace=0.04,
            width_ratios=[3, 1],
            height_ratios=[1, 3],
        )

        self.fig = fig
        self.filename = filename
        self.main_axis = plt.subplot(gs[1, 0])
        self.axis_x = plt.subplot(gs[0, 0])
        self.axis_y = plt.subplot(gs[1, 1])
        self.cax = fig.add_axes([0.15, 0.945, 0.603, 0.020])

    def plot_main(self, image, params, title=None, beam_position=None):
        ny, nx = image.shape

        intensity_max = image.max()
        intensity_min = image.min()

        nx_half = int(nx / 2)
        ny_half = int(ny / 2)
        width_x = int(nx / 2)
        width_y = int(ny / 2)

        self.axis_x.set_xlim(0, nx)
        self.axis_y.set_ylim(0, ny)

        self.main_axis.set_xlabel(r"Pixel Index X")
        self.main_axis.set_ylabel(r"Pixel Index Y")

        if params.projection.color_cutoff:
            vmax = float(params.projection.color_cutoff)
        else:
            temp_image = remove_pixels_by_intensity(image, percent=0.0045)
            vmax = temp_image.max()
        img = self.main_axis.imshow(
            image,
            cmap="jet",
            aspect="auto",
            origin="lower",
            rasterized=True,
            interpolation="none",
            vmax=vmax,
        )

        plt.colorbar(img, self.cax, orientation="horizontal")

        self.main_axis.set_xlim(nx_half - width_x, nx_half + width_x)
        self.main_axis.set_ylim(ny_half - width_y, ny_half + width_y)

        if len(beam_position) != 2:
            raise ValueError(f"Incorrect beam position: {beam_position}")

        # Mark the beam position with a dot
        beam_x, beam_y = beam_position

        self.main_axis.plot([beam_x], [beam_y], marker="o", ms=1.0, c="white", lw=0)

        # Plot circles around the beam center
        diagonal_width = np.sqrt(width_x**2 + width_y**2)
        radial_distance = diagonal_width / 10
        radii = np.arange(0, 2 * diagonal_width, radial_distance)
        for radius in radii:
            c = Circle(
                (beam_x, beam_y),
                radius,
                facecolor="none",
                edgecolor="white",
                lw=0.7,
                ls=(1, (2, 2)),
            )
            self.main_axis.add_patch(c)

        self.main_axis.axvline(beam_x, lw=0.7, ls=(1, (2, 2)), c="white")
        self.main_axis.axhline(beam_y, lw=0.7, ls=(1, (2, 2)), c="white")

        self.axis_x.tick_params(labelbottom=False)
        self.axis_y.tick_params(labelleft=False)

        label = "beam (x, y):"
        self.axis_x.text(
            1.02,
            1.00,
            label,
            va="top",
            fontsize=8,
            ha="left",
            transform=self.axis_x.transAxes,
        )

        beam_position_str = f"({beam_x:.0f}, {beam_y:.0f})"
        self.axis_x.text(
            1.02,
            0.85,
            beam_position_str,
            va="top",
            fontsize=8,
            ha="left",
            transform=self.axis_x.transAxes,
        )

        label = "dimensions (nx, ny):"
        self.axis_x.text(
            1.02,
            0.70,
            label,
            va="top",
            fontsize=8,
            ha="left",
            transform=self.axis_x.transAxes,
        )

        dimensions = f"({nx:.0f}, {ny:.0f})"
        self.axis_x.text(
            1.02,
            0.55,
            dimensions,
            va="top",
            fontsize=8,
            ha="left",
            transform=self.axis_x.transAxes,
        )

        label = f"Imax: {intensity_max:.0f} "
        self.axis_x.text(
            1.02,
            0.40,
            label,
            va="top",
            fontsize=8,
            ha="left",
            transform=self.axis_x.transAxes,
        )

        label = f"Imin: {intensity_min:.0f} "
        self.axis_x.text(
            1.02,
            0.25,
            label,
            va="top",
            fontsize=8,
            ha="left",
            transform=self.axis_x.transAxes,
        )

        if title:
            self.main_axis.text(
                0.15,
                0.999,
                title,
                va="top",
                ha="left",
                transform=self.fig.transFigure,
                fontsize=8,
            )

        ranges_x = params.projection.exclude_pixel_range_x
        spans_x = convert_range_into_spans(ranges_x)

        ranges_y = params.projection.exclude_pixel_range_y
        spans_y = convert_range_into_spans(ranges_y)

        for span in spans_x:
            self.main_axis.axvspan(span.min, span.max, color="white", alpha=0.3, lw=0)
        for span in spans_y:
            self.main_axis.axhspan(span.min, span.max, color="white", alpha=0.3, lw=0)

    def save_and_close(self):
        plt.savefig(self.filename, dpi=600)
        plt.close(self.fig)
