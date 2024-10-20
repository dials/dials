import matplotlib.pyplot as plt
# from collections import namedtuple
# from matplotlib.ticker import MultipleLocator
from matplotlib import gridspec
import numpy as np
from matplotlib.patches import Circle


class Figure:

    def __init__(self, filename='fig.png'):

        fig = plt.figure(figsize=(5, 5))
        gs = gridspec.GridSpec(2, 2, top=0.90, bottom=0.09, left=0.15,
                               right=0.96, wspace=0.04, hspace=0.04,
                               width_ratios=[3, 1], height_ratios=[1, 3])

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

        # mloc = MultipleLocator(0.1)
        # self.axis_x.yaxis.set_minor_locator(mloc)

        # self.axis_x.set_yticks([0, 0.5, 1.0])
        # self.axis_y.set_xticks([0, 0.5, 1.0])

        # self.axis_x.set_ylim(-0.20, 1.1)
        # self.axis_y.set_xlim(-0.20, 1.1)

        self.main_axis.set_xlabel(r"Pixel Index X")
        self.main_axis.set_ylabel(r"Pixel Index Y")

        if params.projection.color_cutoff:
            vmax = float(params.projection.color_cutoff)
        else:
            vmax = image.max()
        img = self.main_axis.imshow(image, cmap="jet", aspect="auto",
                                    origin="lower", rasterized=True,
                                    interpolation="none",
                                    vmax=vmax)

        plt.colorbar(img, self.cax, orientation='horizontal')

        self.main_axis.set_xlim(nx_half - width_x, nx_half + width_x)
        self.main_axis.set_ylim(ny_half - width_y, ny_half + width_y)

        if len(beam_position) != 2:
            raise ValueError(f"Incorrect beam position: {beam_position}")

        # Mark the beam position with a dot
        beam_x, beam_y = beam_position

        self.main_axis.plot([beam_x], [beam_y], marker="o", ms=1.0,
                            c="white", lw=0)

        # Plot circles around the beam center
        diagonal_width = np.sqrt(width_x**2 + width_y**2)
        radial_distance = diagonal_width / 10
        radii = np.arange(0, 2*diagonal_width, radial_distance)
        for radius in radii:
            c = Circle((beam_x, beam_y), radius, facecolor="none",
                       edgecolor="white", lw=0.7, ls=(1, (2, 2)))
            self.main_axis.add_patch(c)

        self.main_axis.axvline(beam_x, lw=0.7, ls=(1, (2, 2)), c='white')
        self.main_axis.axhline(beam_y, lw=0.7, ls=(1, (2, 2)), c='white')

        self.axis_x.tick_params(labelbottom=False)
        self.axis_y.tick_params(labelleft=False)
        self.axis_y.tick_params(axis='x', colors='blue')
        self.axis_x.tick_params(axis='y', colors='blue')

        label = "beam (x, y):"
        self.axis_x.text(1.02, 1.00, label, va='top', fontsize=8,
                         ha='left', transform=self.axis_x.transAxes)

        beam_position_str = f"({beam_x:.0f}, {beam_y:.0f})"
        self.axis_x.text(1.02, 0.85, beam_position_str, va='top', fontsize=8,
                         ha='left', transform=self.axis_x.transAxes)

        label = "dimensions (nx, ny):"
        self.axis_x.text(1.02, 0.70, label, va='top', fontsize=8,
                         ha='left', transform=self.axis_x.transAxes)

        dimensions = f"({nx:.0f}, {ny:.0f})"
        self.axis_x.text(1.02, 0.55, dimensions, va='top', fontsize=8,
                         ha='left', transform=self.axis_x.transAxes)

        label = f"Imax: {intensity_max:.0f} "
        self.axis_x.text(1.02, 0.40, label, va='top', fontsize=8,
                         ha='left', transform=self.axis_x.transAxes)

        label = f"Imin: {intensity_min:.0f} "
        self.axis_x.text(1.02, 0.25, label, va='top', fontsize=8,
                         ha='left', transform=self.axis_x.transAxes)

        if title:
            self.main_axis.text(0.15, 0.999, title, va='top', ha='left',
                                transform=self.fig.transFigure, fontsize=8)

    def save_and_close(self):

        plt.savefig(self.filename, dpi=600)
        plt.close(self.fig)


#    if params.span_xy is not None:
#        if len(params.span_xy) != 4:
#            raise ValueError("span_xy should be a tuple of four numbers")
#        i1, i2, j1, j2 = params.span_xy
#        ax_x.axvspan(i1, i2, color="#BEBEBE", alpha=0.5)
#        ax_y.axhspan(j1, j2, color="#BEBEBE", alpha=0.5)

    # Plot the diffraction image
#
#    bx, by = params.beam_position
#    ax_x.text(0.0, 1.05, f"({bx:.0f}, {by:.0f})",
#              transform=ax_x.transAxes)
#    ax_x.text(0.4, 1.05, f"max: {params.image.max():.2f}",
#              transform=ax_x.transAxes)
#    ax_x.text(0.9, 1.05, f"avg: {params.image.mean():.2f}",
#              transform=ax_x.transAxes)
#
#    # Plot projected profiles
#    for px in params.profiles_x:
#        ax_x.plot(px.x, px.y, c=px.c, lw=px.lw, ms=px.ms, marker=px.marker)
#    for py in params.profiles_y:
#        ax_y.plot(py.x, py.y, c=py.c, lw=py.lw, ms=py.ms, marker=py.marker)
#
#    decorate_plot(params.image, ax, ax_x, ax_y, params.beam_position)
#
#    print(f"Filename '{params.filename}'")
#    plt.savefig(params.filename, dpi=400)
#    plt.close(fig)
