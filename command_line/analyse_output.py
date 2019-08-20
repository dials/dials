#!/usr/bin/env python
#
# analyse_output.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import copy
import errno
import os
import math

import matplotlib

import libtbx.phil
from dials.array_family import flex

# Offline backend
matplotlib.use("Agg")

from matplotlib import pyplot

RAD2DEG = 180 / math.pi

help_message = """

Generate a number of analysis plots from input integrated or indexed reflections.

Examples::

  dials.analyse_output indexed.refl

  dials.analyse_output refined.refl

  dials.analyse_output integrated.refl

"""

# Create the phil parameters
phil_scope = libtbx.phil.parse(
    """
  output {
    directory = .
      .type = str
      .help = "The directory to store the results"
  }
  grid_size = Auto
    .type = ints(size=2)
  pixels_per_bin = 10
    .type = int(value_min=1)

  centroid_diff_max = None
    .help = "Magnitude in pixels of shifts mapped to the extreme colours"
            "in the heatmap plots centroid_diff_x and centroid_diff_y"
    .type = float
    .expert_level = 1
"""
)


def ensure_directory(path):
    """ Make the directory if not already there. """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def ensure_required(rlist, required):
    """ Check which keys aren't present. """
    not_present = []
    for k in required:
        if k not in rlist:
            not_present.append(k)
    if len(not_present) != 0:
        print(" Skipping: following required fields not present:")
        for k in not_present:
            print("  %s" % k)
        return False
    return True


def determine_grid_size(rlist, grid_size=None):
    from libtbx import Auto

    panel_ids = rlist["panel"]
    n_panels = flex.max(panel_ids) + 1
    if grid_size is not None and grid_size is not Auto:
        assert (grid_size[0] * grid_size[1]) >= n_panels, n_panels
        return grid_size
    n_cols = int(math.floor(math.sqrt(n_panels)))
    n_rows = int(math.ceil(n_panels / n_cols))
    return n_cols, n_rows


class per_panel_plot(object):

    title = None
    filename = None
    cbar_ylabel = None
    xlabel = "x"
    ylabel = "y"

    def __init__(self, rlist, directory, grid_size=None, pixels_per_bin=10):
        min_x, max_x, min_y, max_y = self.get_min_max_xy(rlist)
        panel_ids = rlist["panel"]
        crystal_ids = rlist["id"]
        n_crystals = flex.max(crystal_ids) + 1
        n_panels = flex.max(panel_ids) + 1

        n_cols, n_rows = determine_grid_size(rlist, grid_size=grid_size)

        for i_crystal in range(n_crystals):
            crystal_sel = crystal_ids == i_crystal
            fig, axes = pyplot.subplots(n_rows, n_cols, squeeze=False)

            self.gridsize = tuple(
                int(math.ceil(i))
                for i in (max_x / pixels_per_bin, max_y / pixels_per_bin)
            )

            clim = (1e8, 1e-8)

            plots = []

            i_panel = 0
            for i_row in range(n_rows):
                for i_col in range(n_cols):

                    panel_sel = panel_ids == i_panel
                    sel = panel_sel & crystal_sel
                    i_panel += 1

                    if n_panels > 1:
                        axes[i_row][i_col].set_title("Panel %d" % i_panel)
                        axes[i_row][i_col].set_title("Panel %d" % i_panel)

                    if (i_row + 1) == n_rows:
                        axes[i_row][i_col].set_xlabel(self.xlabel)
                    else:
                        pyplot.setp(axes[i_row][i_col].get_xticklabels(), visible=False)

                    if i_col == 0:
                        axes[i_row][i_col].set_ylabel(self.ylabel)
                    else:
                        pyplot.setp(axes[i_row][i_col].get_yticklabels(), visible=False)

                    if sel.count(True) > 0:
                        rlist_sel = rlist.select(sel)
                        if len(rlist_sel) <= 1:
                            ax = pyplot.scatter([], [])  # create empty plot
                        else:
                            ax = self.plot_one_panel(axes[i_row][i_col], rlist_sel)
                            clim = (
                                min(clim[0], ax.get_clim()[0]),
                                max(clim[1], ax.get_clim()[1]),
                            )
                        plots.append(ax)

                    axes[i_row][i_col].set_xlim(min_x, max_x)
                    axes[i_row][i_col].set_ylim(min_y, max_y)
                    axes[i_row][i_col].axes.set_aspect("equal")
                    axes[i_row][i_col].invert_yaxis()

            for p in plots:
                p.set_clim(clim)

            default_size = fig.get_size_inches()
            if self.cbar_ylabel is not None and (n_cols, n_rows) == (1, 24):
                fig.set_size_inches(
                    (n_cols * default_size[0], 0.15 * n_rows * default_size[1])
                )
            elif self.cbar_ylabel is not None and (n_cols, n_rows) == (5, 24):
                fig.set_size_inches(
                    (n_cols * default_size[0], 0.5 * n_rows * default_size[1])
                )
            else:
                fig.set_size_inches(
                    (n_cols * default_size[0], n_rows * default_size[1])
                )

            # pyplot.tight_layout()
            if self.cbar_ylabel is not None:
                cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
                for ax in plots:
                    try:
                        cbar = fig.colorbar(ax, cax=cax)
                        cbar.ax.set_ylabel(self.cbar_ylabel, fontsize=n_cols * 10)
                        cbar.ax.tick_params(labelsize=n_cols * 8)
                    except Exception:
                        continue
                    else:
                        break
                if 1 and (n_cols, n_rows) == (1, 24):
                    fig.subplots_adjust(hspace=0.1 / (n_rows), right=0.8)
                elif n_panels > 1:
                    fig.subplots_adjust(hspace=0.1 / n_rows, right=0.8)

            if self.title is not None:
                fig.suptitle(self.title, fontsize=n_cols * 12)
            fig.savefig(os.path.join(directory, self.filename))
            fig.set_size_inches(default_size)
            pyplot.close()

    def get_min_max_xy(self, rlist):
        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()

        min_x = math.floor(min(flex.min(xc), flex.min(xo)))
        min_y = math.floor(min(flex.min(yc), flex.min(yo)))
        max_x = math.ceil(max(flex.max(xc), flex.max(xo)))
        max_y = math.ceil(max(flex.max(yc), flex.max(yo)))
        return min_x, max_x, min_y, max_y

    def plot_one_panel(self, ax, rlist):
        raise NotImplementedError()


class StrongSpotsAnalyser(object):
    """ Analyse a list of strong spots. """

    def __init__(self, directory):
        """ Setup the directory. """

        # Set the directory
        self.directory = os.path.join(directory, "strong")
        ensure_directory(self.directory)

        # Set the required fields
        self.required = ["xyzobs.px.value", "panel"]

    def __call__(self, rlist):
        """ Analyse the strong spots. """
        from dials.util.command_line import Command

        # Check we have the required fields
        print("Analysing strong spots")
        if not ensure_required(rlist, self.required):
            return

        # Remove I_sigma <= 0
        selection = rlist["intensity.sum.variance"] <= 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(" Removing %d reflections with variance <= 0" % selection.count(True))

        if "flags" in rlist:
            # Select only strong reflections
            Command.start(" Selecting only strong reflections")
            mask = rlist.get_flags(rlist.flags.strong)
            if mask.count(True) > 0:
                rlist = rlist.select(mask)
            Command.end(" Selected %d strong reflections" % len(rlist))

        # Look at distribution of spot counts
        self.spot_count_per_image(rlist)
        self.spot_count_per_panel(rlist)

    def spot_count_per_image(self, rlist):
        """ Analyse the spot count per image. """
        x, y, z = rlist["xyzobs.px.value"].parts()
        max_z = int(math.ceil(flex.max(z)))

        ids = rlist["id"]
        spot_count_per_image = []
        for j in range(flex.max(ids) + 1):
            spot_count_per_image.append(flex.int())
            zsel = z.select(ids == j)
            for i in range(max_z):
                sel = (zsel >= i) & (zsel < (i + 1))
                spot_count_per_image[j].append(sel.count(True))

        colours = ["blue", "red", "green", "orange", "purple", "black"] * 10
        assert len(spot_count_per_image) <= colours

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_title("Spot count per image")
        for j in range(len(spot_count_per_image)):
            ax.scatter(
                list(range(len(spot_count_per_image[j]))),
                spot_count_per_image[j],
                s=5,
                color=colours[j],
                marker="o",
                alpha=0.4,
            )
        ax.set_xlabel("Image #")
        ax.set_ylabel("# spots")
        pyplot.savefig(os.path.join(self.directory, "spots_per_image.png"))
        pyplot.close()

    def spot_count_per_panel(self, rlist):
        """ Analyse the spot count per panel. """
        panel = rlist["panel"]
        if flex.max(panel) == 0:
            # only one panel, don't bother generating a plot
            return

        n_panels = int(flex.max(panel))
        spot_count_per_panel = flex.int()
        for i in range(n_panels):
            sel = (panel >= i) & (panel < (i + 1))
            spot_count_per_panel.append(sel.count(True))

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_title("Spot count per panel")
        ax.scatter(
            list(range(len(spot_count_per_panel))),
            spot_count_per_panel,
            s=10,
            color="blue",
            marker="o",
            alpha=0.4,
        )
        ax.set_xlabel("Panel #")
        ax.set_ylabel("# spots")
        pyplot.savefig(os.path.join(self.directory, "spots_per_panel.png"))
        pyplot.close()


class CentroidAnalyser(object):
    """ Analyse the reflection centroids. """

    def __init__(
        self, directory, grid_size=None, pixels_per_bin=10, centroid_diff_max=1.5
    ):
        """ Setup the directory. """
        # Set the directory
        self.directory = os.path.join(directory, "centroid")
        ensure_directory(self.directory)
        self.grid_size = grid_size
        self.pixels_per_bin = pixels_per_bin
        self.centroid_diff_max = centroid_diff_max

        # Set the required fields
        self.required = [
            "intensity.sum.value",
            "intensity.sum.variance",
            "xyzcal.px",
            "xyzobs.px.value",
            "xyzcal.mm",
            "xyzobs.mm.value",
        ]

    def __call__(self, rlist):
        """ Analyse the reflection centroids. """
        from dials.util.command_line import Command

        # Check we have the required fields
        print("Analysing reflection centroids")
        if not ensure_required(rlist, self.required):
            return

        # Remove I_sigma <= 0
        selection = rlist["intensity.sum.variance"] <= 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(" Removing %d reflections with variance <= 0" % selection.count(True))

        # Remove partial reflections as their observed centroids won't be accurate
        if "partiality" in rlist:
            selection = rlist["partiality"] < 0.99
            if selection.count(True) > 0 and selection.count(True) < selection.size():
                rlist.del_selected(selection)
                print(" Removing %d partial reflections" % selection.count(True))

        # Select only integrated reflections
        Command.start(" Selecting only summation-integated reflections")
        mask = rlist.get_flags(rlist.flags.integrated_sum)
        if mask.count(True) > 0:
            threshold = 10
            rlist = rlist.select(mask)
            Command.end(" Selected %d summation-integrated reflections" % len(rlist))
        else:
            # Select only those reflections used in refinement
            threshold = 0
            mask = rlist.get_flags(rlist.flags.used_in_refinement)
            rlist = rlist.select(mask)
            Command.end(" Selected %d refined reflections" % len(rlist))

        # Look at differences in calculated/observed position
        print(" Analysing centroid differences with I/Sigma > %s" % threshold)
        self.centroid_diff_hist(rlist, threshold)
        print(" Analysing centroid differences in x/y with I/Sigma > %s" % threshold)
        self.centroid_diff_xy(rlist, threshold)
        self.centroid_xy_xz_zy_residuals(rlist, threshold)
        print(" Analysing centroid differences in z with I/Sigma > %s" % threshold)
        self.centroid_diff_z(rlist, threshold)
        print(" Analysing centroid differences vs phi with I/Sigma > %s" % threshold)
        self.centroid_mean_diff_vs_phi(rlist, threshold)

    def centroid_diff_hist(self, rlist, threshold):
        """ Analyse the correlations. """
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        rlist = rlist.select(mask)
        assert len(rlist) > 0
        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()
        xd = xo - xc
        yd = yo - yc
        zd = zo - zc
        diff = flex.sqrt(xd * xd + yd * yd + zd * zd)
        fig = pyplot.figure()
        pyplot.title("Difference between observed and calculated")
        pyplot.hist(diff, bins=20)
        pyplot.xlabel("Difference in position (pixels)")
        pyplot.ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "centroid_diff_hist.png"))
        pyplot.close()

    def centroid_diff_xy(self, rlist, threshold):
        """ Look at the centroid difference in x, y """
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        rlist = rlist.select(mask)
        assert len(rlist) > 0

        class diff_x_plot(per_panel_plot):
            def __init__(self, *args, **kwargs):

                self.title = "Difference between observed and calculated in X"
                self.filename = "centroid_diff_x.png"
                self.cbar_ylabel = "Difference in x position (pixels)"
                self.centroid_diff_max = kwargs.pop("centroid_diff_max", None)
                super(diff_x_plot, self).__init__(*args, **kwargs)

            def plot_one_panel(self, ax, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                xd = xo - xc

                if self.centroid_diff_max is None:
                    self.centroid_diff_max = max(abs(xd))

                hex_ax = ax.hexbin(
                    xc.as_numpy_array(),
                    yc.as_numpy_array(),
                    C=xd.as_numpy_array(),
                    gridsize=self.gridsize,
                    vmin=-1.0 * self.centroid_diff_max,
                    vmax=self.centroid_diff_max,
                )
                return hex_ax

        class diff_y_plot(per_panel_plot):
            def __init__(self, *args, **kwargs):

                self.title = "Difference between observed and calculated in Y"
                self.filename = "centroid_diff_y.png"
                self.cbar_ylabel = "Difference in y position (pixels)"
                self.centroid_diff_max = kwargs.pop("centroid_diff_max", None)
                super(diff_y_plot, self).__init__(*args, **kwargs)

            def plot_one_panel(self, ax, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                yd = yo - yc

                if self.centroid_diff_max is None:
                    self.centroid_diff_max = max(abs(yd))

                hex_ax = ax.hexbin(
                    xc.as_numpy_array(),
                    yc.as_numpy_array(),
                    C=yd.as_numpy_array(),
                    gridsize=self.gridsize,
                    vmin=-1.0 * self.centroid_diff_max,
                    vmax=self.centroid_diff_max,
                )
                return hex_ax

        diff_x_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
            centroid_diff_max=self.centroid_diff_max,
        )
        diff_y_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
            centroid_diff_max=self.centroid_diff_max,
        )

    def centroid_diff_z(self, rlist, threshold):
        """ Look at the centroid difference in x, y """
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        rlist = rlist.select(mask)
        assert len(rlist) > 0
        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()
        zd = zo - zc
        fig = pyplot.figure()
        pyplot.title("Difference between observed and calculated in Z")
        cax = pyplot.hexbin(zc, zd, gridsize=100)
        cax.axes.set_xlabel("z (images)")
        cax.axes.set_ylabel("Difference in z position")
        cbar = pyplot.colorbar(cax)
        cbar.ax.set_ylabel("# Reflections")
        fig.savefig(os.path.join(self.directory, "centroid_diff_z.png"))
        pyplot.close()

    def centroid_mean_diff_vs_phi(self, rlist, threshold):
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        rlist = rlist.select(mask)
        assert len(rlist) > 0

        xc, yc, zc = rlist["xyzcal.mm"].parts()
        xo, yo, zo = rlist["xyzobs.mm.value"].parts()

        dx = xc - xo
        dy = yc - yo
        dphi = (zc - zo) * RAD2DEG

        mean_residuals_x = flex.double()
        mean_residuals_y = flex.double()
        mean_residuals_phi = flex.double()
        rmsd_x = flex.double()
        rmsd_y = flex.double()
        rmsd_phi = flex.double()
        phi_obs_deg = RAD2DEG * zo
        phi = []

        for i_phi in range(
            int(math.floor(flex.min(phi_obs_deg))),
            int(math.ceil(flex.max(phi_obs_deg))),
        ):
            sel = (phi_obs_deg >= i_phi) & (phi_obs_deg < (i_phi + 1))
            if sel.count(True) == 0:
                continue
            mean_residuals_x.append(flex.mean(dx.select(sel)))
            mean_residuals_y.append(flex.mean(dy.select(sel)))
            mean_residuals_phi.append(flex.mean(dphi.select(sel)))
            rmsd_x.append(math.sqrt(flex.mean_sq(dx.select(sel))))
            rmsd_y.append(math.sqrt(flex.mean_sq(dy.select(sel))))
            rmsd_phi.append(math.sqrt(flex.mean_sq(dphi.select(sel))))
            phi.append(i_phi)

        fig = pyplot.figure()
        ax = fig.add_subplot(311)
        # fig.subplots_adjust(hspace=0.5)
        pyplot.axhline(0, color="grey")
        ax.scatter(phi, mean_residuals_x)
        ax.set_xlabel("phi (deg)")
        ax.set_ylabel(r"mean $\Delta$ x (mm)")
        ax = fig.add_subplot(312)
        pyplot.axhline(0, color="grey")
        ax.scatter(phi, mean_residuals_y)
        ax.set_xlabel("phi (deg)")
        ax.set_ylabel(r"mean $\Delta$ y (mm)")
        ax = fig.add_subplot(313)
        pyplot.axhline(0, color="grey")
        ax.scatter(phi, mean_residuals_phi)
        ax.set_xlabel("phi (deg)")
        ax.set_ylabel(r"mean $\Delta$ phi (deg)")
        pyplot.savefig(os.path.join(self.directory, "centroid_mean_diff_vs_phi.png"))
        pyplot.close()

        fig = pyplot.figure()
        ax = fig.add_subplot(311)
        # fig.subplots_adjust(hspace=0.5)
        pyplot.axhline(flex.mean(rmsd_x), color="grey")
        ax.scatter(phi, rmsd_x)
        ax.set_xlabel("phi (deg)")
        ax.set_ylabel("rmsd x (mm)")
        ax = fig.add_subplot(312)
        pyplot.axhline(flex.mean(rmsd_y), color="grey")
        ax.scatter(phi, rmsd_y)
        ax.set_xlabel("phi (deg)")
        ax.set_ylabel("rmsd y (mm)")
        ax = fig.add_subplot(313)
        pyplot.axhline(flex.mean(rmsd_phi), color="grey")
        ax.scatter(phi, rmsd_phi)
        ax.set_xlabel("phi (deg)")
        ax.set_ylabel("rmsd phi (deg)")
        pyplot.savefig(os.path.join(self.directory, "centroid_rmsd_vs_phi.png"))
        pyplot.close()

    def centroid_xy_xz_zy_residuals(self, rlist, threshold):
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        rlist = rlist.select(mask)
        assert len(rlist) > 0

        class residuals_xy_plot(per_panel_plot):

            title = "Centroid residuals in X and Y"
            filename = "centroid_xy_residuals.png"
            cbar_ylabel = None
            xlabel = "X (pixels)"
            ylabel = "Y (pixels)"

            def plot_one_panel(self, ax, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                dx = xc - xo
                dy = yc - yo

                ax.axhline(0, color="grey")
                ax.axvline(0, color="grey")
                ax_xy = ax.scatter(
                    dx.as_numpy_array(), dy.as_numpy_array(), c="b", alpha=0.3
                )
                ax.set_aspect("equal")
                return ax_xy

            def get_min_max_xy(self, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                dx = xc - xo
                dy = yc - yo

                min_x = math.floor(flex.min(dx))
                min_y = math.floor(flex.min(dy))
                max_x = math.ceil(flex.max(dx))
                max_y = math.ceil(flex.max(dy))
                return min_x, max_x, min_y, max_y

        class residuals_zy_plot(per_panel_plot):

            title = "Centroid residuals in Z and Y"
            filename = "centroid_zy_residuals.png"
            cbar_ylabel = None
            xlabel = "Z (images)"
            ylabel = "Y (pixels)"

            def plot_one_panel(self, ax, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                dy = yc - yo
                dz = zc - zo

                ax.axhline(0, color="grey")
                ax.axvline(0, color="grey")
                ax_zy = ax.scatter(
                    dz.as_numpy_array(), dy.as_numpy_array(), c="b", alpha=0.3
                )
                ax.set_aspect("equal")

                return ax_zy

            def get_min_max_xy(self, rlist):
                _, yc, zc = rlist["xyzcal.px"].parts()
                _, yo, zo = rlist["xyzobs.px.value"].parts()
                dy = yc - yo
                dz = zc - zo

                min_x = math.floor(flex.min(dz))
                min_y = math.floor(flex.min(dy))
                max_x = math.ceil(flex.max(dz))
                max_y = math.ceil(flex.max(dy))
                return min_x, max_x, min_y, max_y

        class residuals_xz_plot(per_panel_plot):

            title = "Centroid residuals in X and Z"
            filename = "centroid_xz_residuals.png"
            cbar_ylabel = None
            xlabel = "X (pixels)"
            ylabel = "Z (images)"

            def plot_one_panel(self, ax, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                dx = xc - xo
                dz = zc - zo

                ax.axhline(0, color="grey")
                ax.axvline(0, color="grey")
                ax_xz = ax.scatter(
                    dx.as_numpy_array(), dz.as_numpy_array(), c="b", alpha=0.3
                )
                ax.set_aspect("equal")

                return ax_xz

            def get_min_max_xy(self, rlist):
                xc, yc, zc = rlist["xyzcal.px"].parts()
                xo, yo, zo = rlist["xyzobs.px.value"].parts()
                dx = xc - xo
                dz = zc - zo

                min_x = math.floor(flex.min(dx))
                min_y = math.floor(flex.min(dz))
                max_x = math.ceil(flex.max(dx))
                max_y = math.ceil(flex.max(dz))
                return min_x, max_x, min_y, max_y

        residuals_xy_plot(rlist, self.directory, grid_size=self.grid_size)
        residuals_zy_plot(rlist, self.directory, grid_size=self.grid_size)
        residuals_xz_plot(rlist, self.directory, grid_size=self.grid_size)


class BackgroundAnalyser(object):
    """ Analyse the background. """

    def __init__(self, directory, grid_size=None, pixels_per_bin=10):
        """ Setup the directory. """
        # Set the directory
        self.directory = os.path.join(directory, "background")
        ensure_directory(self.directory)
        self.grid_size = grid_size
        self.pixels_per_bin = pixels_per_bin

        # Set the required fields
        self.required = [
            "background.mse",
            "background.mean",
            "intensity.sum.value",
            "intensity.sum.variance",
            "xyzcal.px",
        ]

    def __call__(self, rlist):
        """ Analyse the relfection background. """
        from dials.util.command_line import Command

        # Check we have the required fields
        print("Analysing reflection backgrounds")
        if not ensure_required(rlist, self.required):
            return

        selection = rlist["intensity.sum.variance"] <= 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(" Removing %d reflections with variance <= 0" % selection.count(True))

        selection = rlist["background.mse"] < 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(
                " Removing %d reflections with negative background model RMSD"
                % selection.count(True)
            )

        selection = rlist["background.mean"] <= 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(
                " Removing %d reflections with mean background <= 0"
                % selection.count(True)
            )

        # Select only integrated reflections
        Command.start(" Selecting only integated reflections")
        mask = rlist.get_flags(rlist.flags.integrated)
        if mask.count(True) == 0:
            return

        rlist = rlist.select(mask)
        Command.end(" Selected %d integrated reflections" % len(rlist))

        # Look at distribution of I/Sigma
        print(" Analysing distribution of background mean")
        self.mean_hist(rlist)
        print(" Analysing distribution of background mean vs XY")
        self.mean_vs_xy(rlist)
        print(" Analysing distribution of background mean vs z")
        self.mean_vs_z(rlist)
        print(" Analysing distribution of background mean vs I/Sigma")
        self.mean_vs_ios(rlist)
        print(" Analysing distribution of background CVRMSD")
        self.rmsd_hist(rlist)
        print(" Analysing distribution of background CVRMSD vs XY")
        self.rmsd_vs_xy(rlist)
        print(" Analysing distribution of background CVRMSD vs z")
        self.rmsd_vs_z(rlist)
        print(" Analysing distribution of background CVRMSD vs I/Sigma")
        self.rmsd_vs_ios(rlist)

    def mean_hist(self, rlist):
        """ Analyse the background RMSD. """
        MEAN = rlist["background.mean"]
        fig = pyplot.figure()
        pyplot.title("Background Model mean histogram")
        pyplot.hist(MEAN, bins=20)
        pyplot.xlabel("mean")
        pyplot.ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "background_model_mean_hist"))
        pyplot.close()

    def mean_vs_xy(self, rlist):
        """ Plot I/Sigma vs X/Y """

        class mean_vs_xy_plot(per_panel_plot):

            title = "Distribution of Background Model mean vs X/Y"
            filename = "background_model_mean_vs_xy.png"
            cbar_ylabel = "Background Model mean"

            def plot_one_panel(self, ax, rlist):
                MEAN = rlist["background.mean"]
                x, y, z = rlist["xyzcal.px"].parts()

                hex_ax = ax.hexbin(
                    x.as_numpy_array(),
                    y.as_numpy_array(),
                    C=MEAN.as_numpy_array(),
                    gridsize=self.gridsize,
                    vmin=0,
                    vmax=1,
                )
                return hex_ax

        mean_vs_xy_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
        )

    def mean_vs_z(self, rlist):
        """ Plot I/Sigma vs Z. """
        MEAN = rlist["background.mean"]
        x, y, z = rlist["xyzcal.px"].parts()
        fig = pyplot.figure()
        pyplot.title("Distribution of Background Model mean vs Z")
        cax = pyplot.hexbin(z, MEAN, gridsize=100)
        cax.axes.set_xlabel("z (images)")
        cax.axes.set_ylabel("Background Model mean")
        cbar = pyplot.colorbar(cax)
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "background_model_mean_vs_z.png"))
        pyplot.close()

    def mean_vs_ios(self, rlist):
        """ Analyse the correlations. """
        MEAN = rlist["background.mean"]
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > 0.1
        I_over_S = I_over_S.select(mask)
        MEAN = MEAN.select(mask)
        fig = pyplot.figure()
        pyplot.title("Background Model mean vs Log I/Sigma")
        cax = pyplot.hexbin(flex.log(I_over_S), MEAN, gridsize=100)
        cbar = pyplot.colorbar(cax)
        cax.axes.set_xlabel("Log I/Sigma")
        cax.axes.set_ylabel("Background Model mean")
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "background_model_mean_vs_ios.png"))
        pyplot.close()

    def rmsd_hist(self, rlist):
        """ Analyse the background RMSD. """
        RMSD = flex.sqrt(rlist["background.mse"])
        MEAN = rlist["background.mean"]
        RMSD = RMSD / MEAN
        fig = pyplot.figure()
        pyplot.title("Background Model mean histogram")
        pyplot.hist(RMSD, bins=20)
        pyplot.xlabel("mean")
        pyplot.ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "background_model_cvrmsd_hist"))
        pyplot.close()

    def rmsd_vs_xy(self, rlist):
        """ Plot I/Sigma vs X/Y """

        class rmsd_vs_xy_plot(per_panel_plot):

            title = "Distribution of Background Model CVRMSD vs X/Y"
            filename = "background_model_cvrmsd_vs_xy.png"
            cbar_ylabel = "Background Model CVRMSD"

            def plot_one_panel(self, ax, rlist):
                RMSD = flex.sqrt(rlist["background.mse"])
                MEAN = rlist["background.mean"]
                RMSD = RMSD / MEAN
                x, y, z = rlist["xyzcal.px"].parts()

                hex_ax = ax.hexbin(
                    x.as_numpy_array(),
                    y.as_numpy_array(),
                    C=RMSD.as_numpy_array(),
                    gridsize=self.gridsize,
                    vmin=0,
                    vmax=1,
                )
                return hex_ax

        rmsd_vs_xy_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
        )

    def rmsd_vs_z(self, rlist):
        """ Plot I/Sigma vs Z. """
        RMSD = flex.sqrt(rlist["background.mse"])
        MEAN = rlist["background.mean"]
        RMSD = RMSD / MEAN
        x, y, z = rlist["xyzcal.px"].parts()
        fig = pyplot.figure()
        pyplot.title("Distribution of Background Model CVRMSD vs Z")
        cax = pyplot.hexbin(z, RMSD, gridsize=100)
        cax.axes.set_xlabel("z (images)")
        cax.axes.set_ylabel("Background Model CVRMSD")
        cbar = pyplot.colorbar(cax)
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "background_model_cvrmsd_vs_z.png"))
        pyplot.close()

    def rmsd_vs_ios(self, rlist):
        """ Analyse the correlations. """
        RMSD = flex.sqrt(rlist["background.mse"])
        MEAN = rlist["background.mean"]
        RMSD = RMSD / MEAN
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > 0.1
        I_over_S = I_over_S.select(mask)
        RMSD = RMSD.select(mask)
        fig = pyplot.figure()
        pyplot.title("Background Model CVRMSD vs Log I/Sigma")
        cax = pyplot.hexbin(flex.log(I_over_S), RMSD, gridsize=100)
        cbar = pyplot.colorbar(cax)
        cax.axes.set_xlabel("Log I/Sigma")
        cax.axes.set_ylabel("Background Model CVRMSD")
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "background_model_cvrmsd_vs_ios.png"))
        pyplot.close()


class IntensityAnalyser(object):
    """ Analyse the intensities. """

    def __init__(self, directory, grid_size=None, pixels_per_bin=10):
        """ Set up the directory. """
        # Set the directory
        self.directory = os.path.join(directory, "intensity")
        ensure_directory(self.directory)
        self.grid_size = grid_size
        self.pixels_per_bin = pixels_per_bin

        # Set the required fields
        self.required = ["intensity.sum.value", "intensity.sum.variance", "xyzcal.px"]

    def __call__(self, rlist):
        """ Analyse the reflection centroids. """
        from dials.util.command_line import Command

        # FIXME Do the same and a comparison for intensity.prf

        # Check we have the required fields
        print("Analysing reflection intensities")
        if not ensure_required(rlist, self.required):
            return

        selection = rlist["intensity.sum.variance"] <= 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(" Removing %d reflections with variance <= 0" % selection.count(True))

        selection = rlist["intensity.sum.value"] <= 0
        if selection.count(True) > 0:
            rlist.del_selected(selection)
            print(
                " Removing %d reflections with intensity <= 0" % selection.count(True)
            )

        # Select only integrated reflections
        Command.start(" Selecting only integated reflections")
        mask = rlist.get_flags(rlist.flags.integrated)
        if mask.count(True) == 0:
            return

        rlist = rlist.select(mask)
        Command.end(" Selected %d integrated reflections" % len(rlist))

        # Look at distribution of I/Sigma
        print(" Analysing distribution of I/Sigma")
        self.i_over_s_hist(rlist)
        print(" Analysing distribution of I/Sigma vs xy")
        self.i_over_s_vs_xy(rlist, "sum")
        print(" Analysing distribution of I/Sigma vs xy")
        self.i_over_s_vs_xy(rlist, "prf")
        print(" Analysing distribution of I/Sigma vs z")
        self.i_over_s_vs_z(rlist)
        print(" Analysing number of background pixels used")
        self.num_background_hist(rlist)
        print(" Analysing number of foreground pixels used")
        self.num_foreground_hist(rlist)

    def i_over_s_hist(self, rlist):
        """ Analyse the correlations. """
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        fig = pyplot.figure()
        pyplot.title("Log I/Sigma histogram")
        pyplot.hist(flex.log(I_over_S), bins=20)
        pyplot.xlabel("Log I/Sigma")
        pyplot.ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "ioversigma_hist"))
        pyplot.close()

    def i_over_s_vs_xy(self, rlist, intensity_type):
        """ Plot I/Sigma vs X/Y """

        class i_over_s_vs_xy_plot(per_panel_plot):

            title = "Distribution of I/Sigma vs X/Y"
            filename = "ioversigma_%s_vs_xy.png" % intensity_type
            cbar_ylabel = "Log I/Sigma"

            def plot_one_panel(self, ax, rlist):
                I_sig = flex.sqrt(rlist["intensity.%s.variance" % intensity_type])
                sel = I_sig > 0
                rlist = rlist.select(sel)
                I_sig = I_sig.select(sel)
                I = rlist["intensity.%s.value" % intensity_type]
                I_over_S = I / I_sig
                x, y, z = rlist["xyzcal.px"].parts()

                hex_ax = ax.hexbin(
                    x.as_numpy_array(),
                    y.as_numpy_array(),
                    C=flex.log(I_over_S),
                    gridsize=self.gridsize,
                )
                return hex_ax

        i_over_s_vs_xy_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
        )

    def i_over_s_vs_z(self, rlist):
        """ Plot I/Sigma vs Z. """
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        x, y, z = rlist["xyzcal.px"].parts()
        fig = pyplot.figure()
        pyplot.title("Distribution of I/Sigma vs Z")
        cax = pyplot.hexbin(z, flex.log(I_over_S), gridsize=100)
        cax.axes.set_xlabel("z (images)")
        cax.axes.set_ylabel("Log I/Sigma")
        cbar = pyplot.colorbar(cax)
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "ioversigma_vs_z.png"))
        pyplot.close()

    def num_background_hist(self, rlist):
        """ Analyse the number of background pixels. """
        if "n_background" in rlist:
            N = rlist["n_background"]
            fig = pyplot.figure()
            pyplot.title("Num Background Pixel Histogram")
            pyplot.hist(N, bins=20)
            pyplot.xlabel("Number of pixels")
            pyplot.ylabel("# reflections")
            fig.savefig(os.path.join(self.directory, "n_background_hist.png"))
            pyplot.close()

    def num_foreground_hist(self, rlist):
        """ Analyse the number of foreground pixels. """
        if "n_foreground" in rlist:
            N = rlist["n_foreground"]
            fig = pyplot.figure()
            pyplot.title("Num Foreground Pixel Histogram")
            pyplot.hist(N, bins=20)
            pyplot.xlabel("Number of pixels")
            pyplot.ylabel("# reflections")
            fig.savefig(os.path.join(self.directory, "n_foreground_hist.png"))
            pyplot.close()


class ReferenceProfileAnalyser(object):
    """ Analyse the reference profiles. """

    def __init__(self, directory, grid_size=None, pixels_per_bin=10):
        """ Set up the directory. """
        # Set the directory
        self.directory = os.path.join(directory, "reference")
        ensure_directory(self.directory)
        self.grid_size = grid_size
        self.pixels_per_bin = pixels_per_bin

        # Set the required fields
        self.required = [
            "intensity.prf.value",
            "intensity.prf.variance",
            "xyzcal.px",
            "profile.correlation",
        ]

    def __call__(self, rlist):
        """ Analyse the reference profiles. """
        from dials.util.command_line import Command

        # Check we have the required fields
        print("Analysing reference profiles")
        if not ensure_required(rlist, self.required):
            return

        # Select only integrated reflections
        Command.start(" Selecting only integated reflections")
        mask = rlist.get_flags(rlist.flags.integrated)
        if mask.count(True) == 0:
            return

        rlist = rlist.select(mask)
        Command.end(" Selected %d integrated reflections" % len(rlist))

        # Analyse distribution of reference spots
        print(" Analysing reference profile distribution vs x/y")
        self.reference_xy(rlist)
        print(" Analysing reference profile distribution vs z")
        self.reference_z(rlist)

        # Look at correlations between profiles
        def ideal_correlations(filename, rlist):
            """ Call for reference spots and all reflections. """

            print(" Analysing reflection profile correlations")
            self.ideal_reflection_corr_hist(rlist, filename)

            print(" Analysing reflection profile correlations vs x/y")
            self.ideal_reflection_corr_vs_xy(rlist, filename)

            print(" Analysing reflection profile correlations vs z")
            self.ideal_reflection_corr_vs_z(rlist, filename)

            print(" Analysing reflection profile correlations vs I/Sigma")
            self.ideal_reflection_corr_vs_ios(rlist, filename)

        # Look at correlations between profiles
        def correlations(filename, rlist):
            """ Call for reference spots and all reflections. """

            print(" Analysing reflection profile correlations")
            self.reflection_corr_hist(rlist, filename)

            print(" Analysing reflection profile correlations vs x/y")
            self.reflection_corr_vs_xy(rlist, filename)

            print(" Analysing reflection profile correlations vs z")
            self.reflection_corr_vs_z(rlist, filename)

            print(" Analysing reflection profile correlations vs I/Sigma")
            self.reflection_corr_vs_ios(rlist, filename)

        mask = rlist.get_flags(rlist.flags.reference_spot)
        correlations("reference", rlist.select(mask))
        correlations("reflection", rlist)
        ideal_correlations("reference", rlist.select(mask))
        ideal_correlations("reflection", rlist)

    def reference_xy(self, rlist):
        """ Analyse the distribution of reference profiles. """
        mask = rlist.get_flags(rlist.flags.reference_spot)
        rlist = rlist.select(mask)

        class reference_xy_plot(per_panel_plot):

            title = "Reference profiles binned in X/Y"
            filename = "reference_xy.png"
            cbar_ylabel = "# reflections"

            def plot_one_panel(self, ax, rlist):
                x, y, z = rlist["xyzcal.px"].parts()

                hex_ax = ax.hexbin(
                    x.as_numpy_array(), y.as_numpy_array(), gridsize=self.gridsize
                )
                return hex_ax

        reference_xy_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
        )

    def reference_z(self, rlist):
        """ Analyse the distribution of reference profiles. """
        x, y, z = rlist["xyzcal.px"].parts()
        fig = pyplot.figure()
        pyplot.title("Reference profiles binned in Z")
        pyplot.hist(z, bins=20)
        pyplot.xlabel("z (images)")
        pyplot.ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "reference_z.png"))
        pyplot.close()

    def reflection_corr_hist(self, rlist, filename):
        """ Analyse the correlations. """
        corr = rlist["profile.correlation"]
        fig = pyplot.figure()
        pyplot.title("Reflection correlations histogram")
        pyplot.hist(corr, bins=20)
        pyplot.xlabel("Correlation with reference profile")
        pyplot.ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "%s_corr_hist" % filename))
        pyplot.close()

    def reflection_corr_vs_xy(self, rlist, filename):
        """ Analyse the correlations. """
        tmp_filename = filename

        class corr_vs_xy_plot(per_panel_plot):

            title = "Reflection correlations binned in X/Y"
            filename = "%s_corr_vs_xy.png" % tmp_filename
            cbar_ylabel = "Correlation with reference profile"

            def plot_one_panel(self, ax, rlist):
                corr = rlist["profile.correlation"]
                x, y, z = rlist["xyzcal.px"].parts()

                hex_ax = ax.hexbin(
                    x.as_numpy_array(),
                    y.as_numpy_array(),
                    C=corr.as_numpy_array(),
                    gridsize=self.gridsize,
                    vmin=0,
                    vmax=1,
                )
                return hex_ax

        corr_vs_xy_plot(
            rlist,
            self.directory,
            grid_size=self.grid_size,
            pixels_per_bin=self.pixels_per_bin,
        )

    def reflection_corr_vs_z(self, rlist, filename):
        """ Analyse the correlations. """
        corr = rlist["profile.correlation"]
        x, y, z = rlist["xyzcal.px"].parts()
        fig = pyplot.figure()
        pyplot.title("Reflection correlations vs Z")
        cax = pyplot.hexbin(z, corr, gridsize=100)
        cbar = pyplot.colorbar(cax)
        cax.axes.set_xlabel("z (images)")
        cax.axes.set_ylabel("Correlation with reference profile")
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "%s_corr_vs_z.png" % filename))
        pyplot.close()

    def reflection_corr_vs_ios(self, rlist, filename):
        """ Analyse the correlations. """
        corr = rlist["profile.correlation"]
        I = rlist["intensity.prf.value"]
        I_sig = flex.sqrt(rlist["intensity.prf.variance"])
        mask = I_sig > 0
        I = I.select(mask)
        I_sig = I_sig.select(mask)
        corr = corr.select(mask)
        I_over_S = I / I_sig
        mask = I_over_S > 0.1
        I_over_S = I_over_S.select(mask)
        corr = corr.select(mask)
        fig = pyplot.figure()
        pyplot.title("Reflection correlations vs Log I/Sigma")
        cax = pyplot.hexbin(flex.log(I_over_S), corr, gridsize=100)
        cbar = pyplot.colorbar(cax)
        cax.axes.set_xlabel("Log I/Sigma")
        cax.axes.set_ylabel("Correlation with reference profile")
        cbar.ax.set_ylabel("# reflections")
        fig.savefig(os.path.join(self.directory, "%s_corr_vs_ios.png" % filename))
        pyplot.close()

    def ideal_reflection_corr_hist(self, rlist, filename):
        """ Analyse the correlations. """
        if "correlation.ideal.profile" in rlist:
            corr = rlist["correlation.ideal.profile"]
            pyplot.title("Reflection correlations histogram")
            pyplot.hist(corr, bins=20)
            pyplot.xlabel("Correlation with reference profile")
            pyplot.ylabel("# reflections")
            pyplot.savefig(
                os.path.join(self.directory, "ideal_%s_corr_hist" % filename)
            )
            pyplot.close()

    def ideal_reflection_corr_vs_xy(self, rlist, filename):
        """ Analyse the correlations. """
        if "correlation.ideal.profile" in rlist:
            corr = rlist["correlation.ideal.profile"]
            x, y, z = rlist["xyzcal.px"].parts()
            pyplot.title("Reflection correlations binned in X/Y")
            cax = pyplot.hexbin(x, y, C=corr, gridsize=100, vmin=0.0, vmax=1.0)
            cbar = pyplot.colorbar(cax)
            pyplot.xlabel("x (pixels)")
            pyplot.ylabel("y (pixels)")
            cbar.ax.set_ylabel("Correlation with reference profile")
            pyplot.savefig(
                os.path.join(self.directory, "ideal_%s_corr_vs_xy.png" % filename)
            )
            pyplot.close()

    def ideal_reflection_corr_vs_z(self, rlist, filename):
        """ Analyse the correlations. """
        if "correlation.ideal.profile" in rlist:
            corr = rlist["correlation.ideal.profile"]
            x, y, z = rlist["xyzcal.px"].parts()
            pyplot.title("Reflection correlations vs Z")
            cax = pyplot.hexbin(z, corr, gridsize=100)
            cbar = pyplot.colorbar(cax)
            pyplot.xlabel("z (images)")
            pyplot.ylabel("Correlation with reference profile")
            cbar.ax.set_ylabel("# reflections")
            pyplot.savefig(
                os.path.join(self.directory, "ideal_%s_corr_vs_z.png" % filename)
            )
            pyplot.close()

    def ideal_reflection_corr_vs_ios(self, rlist, filename):
        """ Analyse the correlations. """
        if "correlation.ideal.profile" in rlist:
            corr = rlist["correlation.ideal.profile"]
            I = rlist["intensity.prf.value"]
            I_sig = flex.sqrt(rlist["intensity.prf.variance"])
            mask = I_sig > 0
            I = I.select(mask)
            I_sig = I_sig.select(mask)
            corr = corr.select(mask)
            I_over_S = I / I_sig
            mask = I_over_S > 0.1
            I_over_S = I_over_S.select(mask)
            corr = corr.select(mask)
            pyplot.title("Reflection correlations vs Log I/Sigma")
            cax = pyplot.hexbin(flex.log(I_over_S), corr, gridsize=100)
            cbar = pyplot.colorbar(cax)
            pyplot.xlabel("Log I/Sigma")
            pyplot.ylabel("Correlation with reference profile")
            cbar.ax.set_ylabel("# reflections")
            pyplot.savefig(
                os.path.join(self.directory, "ideal_%s_corr_vs_ios.png" % filename)
            )
            pyplot.close()


def analyse(rlist, directory, grid_size=None, pixels_per_bin=10, centroid_diff_max=1.5):
    """ Setup the analysers. """
    directory = os.path.join(directory, "analysis")
    analysers = [
        StrongSpotsAnalyser(directory),
        CentroidAnalyser(
            directory,
            grid_size=grid_size,
            pixels_per_bin=pixels_per_bin,
            centroid_diff_max=centroid_diff_max,
        ),
        BackgroundAnalyser(
            directory, grid_size=grid_size, pixels_per_bin=pixels_per_bin
        ),
        IntensityAnalyser(
            directory, grid_size=grid_size, pixels_per_bin=pixels_per_bin
        ),
        ReferenceProfileAnalyser(
            directory, grid_size=grid_size, pixels_per_bin=pixels_per_bin
        ),
    ]

    """ Do all the analysis. """
    for a in analysers:
        a(copy.deepcopy(rlist))


def run():
    from dials.util.options import OptionParser

    # Create the parser
    usage = "usage: dials.analyse_output [options] observations.refl"
    parser = OptionParser(
        usage=usage, phil=phil_scope, read_reflections=True, epilog=help_message
    )

    parser.add_option(
        "--xkcd",
        action="store_true",
        dest="xkcd",
        default=False,
        help="Special drawing mode",
    )

    # Parse the command line arguments
    params, options = parser.parse_args(show_diff_phil=True)

    if options.xkcd:
        pyplot.xkcd()

    # Show the help
    if len(params.input.reflections) != 1:
        parser.print_help()
        exit(0)

    # Analyse the reflections
    analyse(
        params.input.reflections[0].data,
        params.output.directory,
        grid_size=params.grid_size,
        pixels_per_bin=params.pixels_per_bin,
        centroid_diff_max=params.centroid_diff_max,
    )


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        run()
    except Exception as e:
        halraiser(e)
