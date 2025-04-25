from __future__ import annotations

import copy
import itertools
import json
import math

import numpy as np

import libtbx.phil
from cctbx import uctbx

import dials.util.log
from dials.algorithms.merging.merge import MergingStatisticsData
from dials.algorithms.scaling.model.model import plot_scaling_models
from dials.algorithms.scaling.observers import make_merging_stats_plots
from dials.algorithms.scaling.scaling_library import (
    merging_stats_from_scaled_array,
    scaled_data_as_miller_array,
)
from dials.algorithms.scaling.scaling_utilities import DialsMergingStatisticsError
from dials.array_family import flex
from dials.report.plots import i_over_sig_i_vs_i_plot
from dials.util import show_mail_handle_errors
from dials.util.command_line import Command

RAD2DEG = 180 / math.pi

help_message = """

Generates a html report given the output of various DIALS programs
(observations.refl and/or models.expt).

Examples::

  dials.report strong.refl

  dials.report indexed.refl

  dials.report refined.refl

  dials.report integrated.refl

  dials.report refined.expt

  dials.report integrated.refl integrated.expt
"""

# Create the phil parameters
phil_scope = libtbx.phil.parse(
    """
  output {
    html = dials.report.html
      .type = path
      .help = "The name of the output html file"
    json = None
      .type = path
      .help = "The name of the optional json file containing the plot data"
    external_dependencies = *remote local embed
      .type = choice
      .help = "Whether to use remote external dependencies (files relocatable"
              "but requires an internet connection), local (does not require"
              "internet connection but files may not be relocatable) or embed"
              "all external dependencies (inflates the html file size)."
  }
  grid_size = Auto
    .type = ints(size=2)
  pixels_per_bin = 40
    .type = int(value_min=1)

  centroid_diff_max = None
    .help = "Magnitude in pixels of shifts mapped to the extreme colours"
            "in the heatmap plots centroid_diff_x and centroid_diff_y"
    .type = float
    .expert_level = 1

  orientation_decomposition
    .help = "Options determining how the orientation matrix"
            "decomposition is done. The axes about which to decompose"
            "the matrix into three rotations are chosen here, as well"
            "as whether the rotations are relative to the reference"
            "orientation, taken from the static crystal model"
  {
    e1 = 1. 0. 0.
      .type = floats(size = 3)

    e2 = 0. 1. 0.
      .type = floats(size = 3)

    e3 = 0. 0. 1.
      .type = floats(size = 3)

    relative_to_static_orientation = True
      .type = bool
  }
"""
)


def ensure_directory(path):
    """Make the directory if not already there."""
    import errno
    from os import makedirs

    try:
        makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def ensure_required(rlist, required):
    """Check which keys aren't present."""
    not_present = []
    for k in required:
        if k not in rlist:
            not_present.append(k)
    if len(not_present) != 0:
        print(" Skipping: following required fields not present:")
        for k in not_present:
            print(f"  {k}")
        return False
    return True


def color_repeats(n=1):
    """Set up a cycle through default Plotly colors, repeating each n times"""

    color_list = (
        ["#1f77b4"] * n  # muted blue
        + ["#ff7f0e"] * n  # safety orange
        + ["#2ca02c"] * n  # cooked asparagus green
        + ["#d62728"] * n  # brick red
        + ["#9467bd"] * n  # muted purple
        + ["#8c564b"] * n  # chestnut brown
        + ["#e377c2"] * n  # raspberry yogurt pink
        + ["#7f7f7f"] * n  # middle gray
        + ["#bcbd22"] * n  # curry yellow-green
        + ["#17becf"] * n  # blue-teal
    )
    return itertools.cycle(color_list)


class ScanVaryingCrystalAnalyser:
    """Analyse a scan-varying crystal."""

    def __init__(self, orientation_decomposition):
        # Decomposition axes
        self._e1 = orientation_decomposition.e1
        self._e2 = orientation_decomposition.e2
        self._e3 = orientation_decomposition.e3
        self._relative_to_static_orientation = (
            orientation_decomposition.relative_to_static_orientation
        )
        self._debug = False

    def __call__(self, experiments):
        """Analyse the strong spots."""
        # Check we have the required fields
        print("Analysing scan-varying crystal model")

        d = {}

        if experiments is not None and len(experiments):
            d.update(self.plot_cell(experiments))
            d.update(self.plot_orientation(experiments))
        return {"scan_varying": d}

    def plot_cell(self, experiments):
        """Analyse the scan-varying cell parameters."""

        # cell plot
        dat = []
        for iexp, exp in enumerate(experiments):
            crystal = exp.crystal
            scan = exp.scan

            if crystal is None:
                print("Ignoring absent crystal model")
                continue

            if crystal.num_scan_points == 0:
                print("Ignoring scan-static crystal")
                continue

            scan_pts = list(range(crystal.num_scan_points))
            cells = [crystal.get_unit_cell_at_scan_point(t) for t in scan_pts]
            cell_params = [e.parameters() for e in cells]
            a, b, c, aa, bb, cc = zip(*cell_params)
            aa = [round(i, ndigits=6) for i in aa]
            bb = [round(i, ndigits=6) for i in bb]
            cc = [round(i, ndigits=6) for i in cc]
            start, stop = scan.get_array_range()
            phi = [scan.get_angle_from_array_index(t) for t in range(start, stop + 1)]
            vol = [e.volume() for e in cells]
            cell_dat = {
                "phi": phi,
                "a": a,
                "b": b,
                "c": c,
                "alpha": aa,
                "beta": bb,
                "gamma": cc,
                "volume": vol,
            }
            if self._debug:
                print(f"Crystal in Experiment {iexp}")
                print("Phi\ta\tb\tc\talpha\tbeta\tgamma\tVolume")
                msg = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}"
                line_dat = zip(phi, a, b, c, aa, bb, cc, vol)
                for line in line_dat:
                    print(msg.format(*line))
            dat.append(cell_dat)

        d = {
            "scan_varying_cell": {
                "data": [],
                "layout": {
                    "title": "Scan-varying cell parameters",
                    "xaxis3": {
                        "domain": [0, 1],
                        "anchor": "y7",
                        "title": "rotation angle (°)",
                    },
                    "xaxis2": {"domain": [0.53, 1], "anchor": "y6"},
                    "xaxis": {"domain": [0, 0.47], "anchor": "y3"},
                    "yaxis7": {"domain": [0.0, 0.2], "anchor": "x3", "nticks": 5},
                    "yaxis6": {"domain": [0.3, 0.5], "anchor": "x2", "nticks": 5},
                    "yaxis5": {"domain": [0.55, 0.75], "anchor": "x2", "nticks": 5},
                    "yaxis4": {"domain": [0.8, 1], "anchor": "x2", "nticks": 5},
                    "yaxis3": {"domain": [0.3, 0.5], "anchor": "x", "nticks": 5},
                    "yaxis2": {"domain": [0.55, 0.75], "anchor": "x", "nticks": 5},
                    "yaxis": {"domain": [0.8, 1], "nticks": 5},
                },
                "help": """\
A plot of the smoothly-varying unit cell parameters determined by scan-varying
refinement. Small variations up to ~1% of the linear cell dimensions are not
uncommon. Especially in cases of radiation damage, an increase of cell volume
with respect to dose is often seen. Multi-turn scans typically show some
oscillatory behaviour. This is not expected to be physically real, but indicates
the refinement algorithm accounting for unmodelled features in the data.
""",
            }
        }

        if len(dat) == 1:
            colors = color_repeats()
        else:
            colors = color_repeats(7)

        for cell_dat in dat:
            d["scan_varying_cell"]["data"].extend(
                [
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["a"],
                        "type": "scatter",
                        "name": "a (Å)",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["b"],
                        "type": "scatter",
                        "name": "b (Å)",
                        "xaxis": "x",
                        "yaxis": "y2",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["c"],
                        "type": "scatter",
                        "name": "c (Å)",
                        "xaxis": "x",
                        "yaxis": "y3",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["alpha"],
                        "type": "scatter",
                        "name": "α (°)",
                        "xaxis": "x2",
                        "yaxis": "y4",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["beta"],
                        "type": "scatter",
                        "name": "β (°)",
                        "xaxis": "x2",
                        "yaxis": "y5",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["gamma"],
                        "type": "scatter",
                        "name": "γ (°)",
                        "xaxis": "x2",
                        "yaxis": "y6",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": cell_dat["phi"],
                        "y": cell_dat["volume"],
                        "type": "scatter",
                        "name": "volume (Å³)",
                        "xaxis": "x3",
                        "yaxis": "y7",
                        "marker": {"color": next(colors)},
                    },
                ]
            )
        if not len(dat):
            return {}
        return d

    def plot_orientation(self, experiments):
        from scitbx import matrix

        from dials.algorithms.refinement.rotation_decomposition import (
            solve_r3_rotation_for_angles_given_axes,
        )

        # orientation plot
        dat = []
        for iexp, exp in enumerate(experiments):
            crystal = exp.crystal
            scan = exp.scan

            if crystal is None:
                print("Ignoring absent crystal model")
                continue

            if crystal.num_scan_points == 0:
                print("Ignoring scan-static crystal")
                continue

            scan_pts = list(range(crystal.num_scan_points))
            start, stop = scan.get_array_range()
            phi = [scan.get_angle_from_array_index(t) for t in range(start, stop + 1)]
            Umats = [matrix.sqr(crystal.get_U_at_scan_point(t)) for t in scan_pts]
            if self._relative_to_static_orientation:
                # factor out static U
                Uinv = matrix.sqr(crystal.get_U()).inverse()
                Umats = [U * Uinv for U in Umats]
            # NB e3 and e1 definitions for the crystal are swapped compared
            # with those used inside the solve_r3_rotation_for_angles_given_axes
            # method
            angles = [
                solve_r3_rotation_for_angles_given_axes(
                    U, self._e3, self._e2, self._e1, deg=True
                )
                for U in Umats
            ]
            phi3, phi2, phi1 = zip(*angles)
            angle_dat = {"phi": phi, "phi3": phi3, "phi2": phi2, "phi1": phi1}
            if self._debug:
                print(f"Crystal in Experiment {iexp}")
                print("Image\tphi3\tphi2\tphi1")
                msg = "{0}\t{1}\t{2}\t{3}"
                line_dat = zip(phi, phi3, phi2, phi1)
                for line in line_dat:
                    print(msg.format(*line))
            dat.append(angle_dat)

        d = {
            "scan_varying_orientation": {
                "data": [],
                "layout": {
                    "title": "Scan-varying orientation parameters",
                    "xaxis3": {
                        "domain": [0, 1],
                        "anchor": "y3",
                        "title": "rotation angle (°)",
                    },
                    "xaxis2": {"domain": [0, 1], "anchor": "y3"},
                    "xaxis": {"domain": [0, 1], "anchor": "y3"},
                    "yaxis3": {"domain": [0, 0.3], "anchor": "x3", "nticks": 5},
                    "yaxis2": {"domain": [0.35, 0.65], "anchor": "x2", "nticks": 5},
                    "yaxis": {"domain": [0.7, 1], "nticks": 5},
                },
                "help": """\
A plot of the smoothly-varying crystal orientation parameters determined by
scan-varying refinement. By default, the crystal orientation matrix is
decomposed into rotations along the orthogonal axes of the laboratory frame.
Alternate axes may be chosen as parameters to the program. Small variations of
a few tenths of a degree are not uncommon. Larger changes may indicate a
problem with the model geometry. Monotonic drift along the axis closest to the
rotation axis may indicate that the goniometer rotation speed is recorded
incorrectly.
""",
            }
        }

        if len(dat) == 1:
            colors = color_repeats()
        else:
            colors = color_repeats(3)
        for ori in dat:
            d["scan_varying_orientation"]["data"].extend(
                [
                    {
                        "x": ori["phi"],
                        "y": ori["phi1"],
                        "type": "scatter",
                        "name": "Φ1 (°)",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": ori["phi"],
                        "y": ori["phi2"],
                        "type": "scatter",
                        "name": "Φ2 (°)",
                        "xaxis": "x2",
                        "yaxis": "y2",
                        "marker": {"color": next(colors)},
                    },
                    {
                        "x": ori["phi"],
                        "y": ori["phi3"],
                        "type": "scatter",
                        "name": "Φ3 (°)",
                        "xaxis": "x3",
                        "yaxis": "y3",
                        "marker": {"color": next(colors)},
                    },
                ]
            )
        if not len(dat):
            return {}
        return d


class StrongSpotsAnalyser:
    """Analyse a list of strong spots."""

    def __init__(self, pixels_per_bin=10):
        self.pixels_per_bin = pixels_per_bin

        # Set the required fields
        self.required = ["xyzobs.px.value", "panel"]

    def __call__(self, rlist):
        """Analyse the strong spots."""
        # Check we have the required fields
        print("Analysing strong spots")
        if not ensure_required(rlist, self.required):
            return {"strong": {}}

        # Remove I_sigma <= 0
        if "intensity.sum.variance" in rlist:
            selection = rlist["intensity.sum.variance"] <= 0
            if selection.count(True) > 0:
                rlist.del_selected(selection)
                print(
                    " Removing %d reflections with variance <= 0"
                    % selection.count(True)
                )

        if "flags" in rlist:
            # Select only strong reflections
            Command.start(" Selecting only strong reflections")
            mask = rlist.get_flags(rlist.flags.strong)
            if mask.count(True) > 0:
                rlist = rlist.select(mask)
            Command.end(f" Selected {len(rlist)} strong reflections")

        x, y, z = rlist["xyzobs.px.value"].parts()
        self.nbinsx, self.nbinsy = tuple(
            int(math.ceil(i))
            for i in (
                flex.max(x) / self.pixels_per_bin,
                flex.max(y) / self.pixels_per_bin,
            )
        )

        # Look at distribution of spot counts
        d = self.spot_count_per_image(rlist)
        # Look at distribution of unindexed spots with detector position
        d.update(self.unindexed_count_xy(rlist))
        # Look at distribution of indexed spots with detector position
        d.update(self.indexed_count_xy(rlist))
        return {"strong": d}

    def spot_count_per_image(self, rlist):
        """Analyse the spot count per image."""
        x, y, z = rlist["xyzobs.px.value"].parts()
        max_z = int(math.ceil(flex.max(z)))

        indexed_sel = rlist.get_flags(rlist.flags.indexed)
        n_indexed = indexed_sel.count(True)

        if "imageset_id" in rlist:
            ids = rlist["imageset_id"]
        else:
            ids = rlist["id"]
        spot_count_per_image = []
        indexed_per_image = []
        for j in range(flex.max(ids) + 1):
            spot_count_per_image.append([])
            ids_sel = ids == j
            zsel = z.select(ids_sel)
            for i in range(max_z):
                sel = (zsel >= i) & (zsel < (i + 1))
                spot_count_per_image[j].append(sel.count(True))
            if n_indexed > 0:
                indexed_per_image.append([])
                zsel = z.select(ids_sel & indexed_sel)
                for i in range(max_z):
                    sel = (zsel >= i) & (zsel < (i + 1))
                    indexed_per_image[j].append(sel.count(True))

        d = {
            "spot_count_per_image": {
                "data": [],
                "layout": {
                    "title": "Spot count per image",
                    "xaxis": {"title": "Image"},
                    "yaxis": {"title": "Spot count", "rangemode": "tozero"},
                },
                "help": """\
A plot of the distribution of total and indexed spot count with respect to image
number. A drop off in spot count towards zero at the end of the scan may be
indicative of radiation damage, whereas a sudden fall followed by a sudden rise
in spot count may suggest that the crystal has moved out of the beam. Systematic
variations in spot count with image number may be a result of unit cell
dimensions, variations in volume of crystal intersecting the beam, or
diffraction anisotropy. Large separation between the total and indexed spot
count shows a significant number of unindexed reflections, which may be the
result of further, unidentified lattices, split reflections, reflections due to
ice rings, or poor spot-finding parameters.
""",
            }
        }
        for j in range(flex.max(ids) + 1):
            d["spot_count_per_image"]["data"].append(
                {
                    "x": list(range(len(spot_count_per_image[j]))),
                    "y": spot_count_per_image[j],
                    "type": "scatter",
                    "name": "#spots",
                    "opacity": 0.4,
                }
            )
            if n_indexed > 0:
                d["spot_count_per_image"]["data"].append(
                    {
                        "x": list(range(len(indexed_per_image[j]))),
                        "y": indexed_per_image[j],
                        "type": "scatter",
                        "name": "#indexed",
                        "opacity": 0.4,
                    }
                )

        if indexed_sel.count(True) > 0 and flex.max(rlist["id"]) > 0:
            # multiple lattices
            ids = rlist["id"]
            indexed_per_lattice_per_image = []
            for j in range(flex.max(ids) + 1):
                indexed_per_lattice_per_image.append([])
                zsel = z.select((ids == j) & indexed_sel)
                for i in range(max_z):
                    sel = (zsel >= i) & (zsel < (i + 1))
                    indexed_per_lattice_per_image[j].append(sel.count(True))

            d["indexed_per_lattice_per_image"] = {
                "data": [],
                "layout": {
                    "title": "Indexed reflections per lattice per image",
                    "xaxis": {"title": "Image"},
                    "yaxis": {"title": "Number of reflections", "rangemode": "tozero"},
                },
            }

            for j in range(flex.max(ids) + 1):
                d["indexed_per_lattice_per_image"]["data"].append(
                    {
                        "x": list(range(len(indexed_per_lattice_per_image[j]))),
                        "y": indexed_per_lattice_per_image[j],
                        "type": "scatter",
                        "name": "Lattice #%i" % (j + 1),
                        "opacity": 0.4,
                    }
                )

        return d

    def unindexed_count_xy(self, rlist):
        """Analyse the spot count in x/y."""
        x, y, z = rlist["xyzobs.px.value"].parts()

        indexed_sel = rlist.get_flags(rlist.flags.indexed)
        if indexed_sel.count(True) == 0 or indexed_sel.count(False) == 0:
            return {}

        x = x.select(~indexed_sel).as_numpy_array()
        y = y.select(~indexed_sel).as_numpy_array()

        H, xedges, yedges = np.histogram2d(x, y, bins=(self.nbinsx, self.nbinsy))

        return {
            "n_unindexed_vs_xy": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": "n_unindexed",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Number of unindexed reflections binned in X/Y",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            }
        }

    def indexed_count_xy(self, rlist):
        """Analyse the indexed spot count in x/y."""
        x, y, z = rlist["xyzobs.px.value"].parts()

        indexed_sel = rlist.get_flags(rlist.flags.indexed)
        if indexed_sel.count(True) == 0 or indexed_sel.count(False) == 0:
            return {}

        x = x.select(indexed_sel).as_numpy_array()
        y = y.select(indexed_sel).as_numpy_array()

        H, xedges, yedges = np.histogram2d(x, y, bins=(self.nbinsx, self.nbinsy))

        return {
            "n_indexed_vs_xy": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": "n_indexed",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Number of indexed reflections binned in X/Y",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            }
        }


class CentroidAnalyser:
    """Analyse the reflection centroids."""

    def __init__(self, grid_size=None, pixels_per_bin=10, centroid_diff_max=1.5):
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
        """Analyse the reflection centroids."""
        # Check we have the required fields
        print("Analysing reflection centroids")
        if not ensure_required(rlist, self.required):
            return {"centroid": {}}

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
        Command.start(" Selecting only summation-integrated reflections")
        mask = rlist.get_flags(rlist.flags.integrated_sum)
        if mask.count(True) > 0:
            threshold = 10
            rlist = rlist.select(mask)
            Command.end(f" Selected {len(rlist)} summation-integrated reflections")
        else:
            # Select only those reflections used in refinement
            threshold = 0
            mask = rlist.get_flags(rlist.flags.used_in_refinement)
            rlist = rlist.select(mask)
            Command.end(f" Selected {len(rlist)} refined reflections")

        # Look at differences in calculated/observed position
        print(f" Analysing centroid differences with I/Sigma > {threshold}")
        d = self.centroid_diff_hist(rlist, threshold)
        print(f" Analysing centroid differences in x/y with I/Sigma > {threshold}")
        d.update(self.centroid_diff_xy(rlist, threshold))
        d.update(self.centroid_xy_xz_zy_residuals(rlist, threshold))
        print(f" Analysing centroid differences in z with I/Sigma > {threshold}")
        d.update(self.centroid_diff_z(rlist, threshold))
        print(f" Analysing centroid differences vs phi with I/Sigma > {threshold}")
        d.update(self.centroid_mean_diff_vs_phi(rlist, threshold))
        return {"centroid": d}

    def centroid_diff_hist(self, rlist, threshold):
        """Analyse the correlations."""
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        if mask.count(True) == 0:
            return {}
        rlist = rlist.select(mask)
        assert len(rlist) > 0
        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()
        xd = xo - xc
        yd = yo - yc
        zd = zo - zc
        diff = flex.sqrt(xd * xd + yd * yd + zd * zd)
        hist = flex.histogram(diff, n_slots=20)

        d = {
            "centroid_difference_histogram": {
                "data": [
                    {
                        "x": list(hist.slot_centers()),
                        "y": list(hist.slots()),
                        "type": "bar",
                        "name": "centroid_differences",
                    }
                ],
                "layout": {
                    "title": "Difference between observed and calculated centroids",
                    "xaxis": {"title": "Difference in position (pixels)"},
                    "yaxis": {"title": "Number of reflections"},
                    "bargap": 0,
                },
            }
        }
        return d

    def centroid_diff_xy(self, rlist, threshold):
        """Look at the centroid difference in x, y"""
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        if mask.count(True) == 0:
            return {}
        rlist = rlist.select(mask)
        assert len(rlist) > 0

        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()
        xd = xo - xc
        yd = yo - yc

        nbinsx, nbinsy = tuple(
            int(math.ceil(i))
            for i in (
                flex.max(xc) / self.pixels_per_bin,
                flex.max(yc) / self.pixels_per_bin,
            )
        )

        xc = xc.as_numpy_array()
        yc = yc.as_numpy_array()
        xd = xd.as_numpy_array()
        yd = yd.as_numpy_array()

        H, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy))
        H1, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy), weights=xd)
        H2, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy), weights=yd)

        nonzeros = np.nonzero(H)
        z1 = np.empty(H.shape)
        z1[:] = np.NAN
        z1[nonzeros] = H1[nonzeros] / H[nonzeros]
        z2 = np.empty(H.shape)
        z2[:] = np.NAN
        z2[nonzeros] = H2[nonzeros] / H[nonzeros]

        return {
            "centroid_differences_x": {
                "data": [
                    {
                        "name": "centroid_differences_x",
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": z1.transpose().tolist(),
                        "type": "heatmap",
                        "colorbar": {
                            "title": "Difference in X position (pixels)",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Difference between observed and calculated centroids in X",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            },
            "centroid_differences_y": {
                "data": [
                    {
                        "name": "centroid_differences_y",
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": z2.transpose().tolist(),
                        "type": "heatmap",
                        "colorbar": {
                            "title": "Difference in Y position (pixels)",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Difference between observed and calculated centroids in Y",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            },
        }

    def centroid_diff_z(self, rlist, threshold):
        """Look at the centroid difference in z"""
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        if mask.count(True) == 0:
            return {}
        rlist = rlist.select(mask)
        assert len(rlist) > 0
        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()
        zd = zo - zc

        if zd.all_approx_equal(zd[0]):
            # probably still images, no z residuals
            return {}

        H, xedges, yedges = np.histogram2d(zc, zd, bins=(100, 100))

        return {
            "centroid_differences_z": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": "centroid_differences_z",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Difference between observed and calculated centroids in Z",
                    "xaxis": {"title": "Z", "showgrid": False},
                    "yaxis": {
                        "title": "Difference in Z position (images)",
                        "showgrid": False,
                    },
                },
            }
        }

    def centroid_mean_diff_vs_phi(self, rlist, threshold):
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        if mask.count(True) == 0:
            return {}
        rlist = rlist.select(mask)
        assert len(rlist) > 0

        xc, yc, zc = rlist["xyzcal.mm"].parts()
        xo, yo, zo = rlist["xyzobs.mm.value"].parts()

        dx = xc - xo
        dy = yc - yo
        dphi = (zc - zo) * RAD2DEG

        if dphi.all_approx_equal(dphi[0]):
            # probably still images, no z residuals
            return {}

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

        d = {
            "centroid_mean_differences_vs_phi": {
                "data": [
                    {
                        "x": list(phi),
                        "y": list(mean_residuals_x),
                        "type": "scatter",
                        "name": "mean_dx",
                    },
                    {
                        "x": list(phi),
                        "y": list(mean_residuals_y),
                        "type": "scatter",
                        "name": "mean_dy",
                        "xaxis": "x2",
                        "yaxis": "y2",
                    },
                    {
                        "x": list(phi),
                        "y": list(mean_residuals_phi),
                        "type": "scatter",
                        "name": "mean_dphi",
                        "xaxis": "x3",
                        "yaxis": "y3",
                    },
                ],
                "layout": {
                    "title": "Difference between observed and calculated centroids vs phi",
                    "yaxis3": {"domain": [0, 0.266]},
                    # 'legend': {'traceorder': 'reversed'},
                    "xaxis3": {"anchor": "y3"},
                    "xaxis2": {"anchor": "y2"},
                    "yaxis2": {"domain": [0.366, 0.633]},
                    "yaxis": {"domain": [0.733, 1]},
                },
            },
            "centroid_rmsd_vs_phi": {
                "data": [
                    {
                        "x": list(phi),
                        "y": list(rmsd_x),
                        "type": "scatter",
                        "name": "rmsd_dx",
                    },
                    {
                        "x": list(phi),
                        "y": list(rmsd_y),
                        "type": "scatter",
                        "name": "rmsd_dy",
                        "xaxis": "x2",
                        "yaxis": "y2",
                    },
                    {
                        "x": list(phi),
                        "y": list(rmsd_phi),
                        "type": "scatter",
                        "name": "rmsd_dphi",
                        "xaxis": "x3",
                        "yaxis": "y3",
                    },
                ],
                "layout": {
                    "title": "RMSD between observed and calculated centroids vs phi",
                    "yaxis3": {"domain": [0, 0.266], "rangemode": "tozero"},
                    # 'legend': {'traceorder': 'reversed'},
                    "xaxis3": {"anchor": "y3"},
                    "xaxis2": {"anchor": "y2"},
                    "yaxis2": {"domain": [0.366, 0.633], "rangemode": "tozero"},
                    "yaxis": {"domain": [0.733, 1], "rangemode": "tozero"},
                },
            },
        }
        return d

    def centroid_xy_xz_zy_residuals(self, rlist, threshold):
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        mask = I_over_S > threshold
        if mask.count(True) == 0:
            return {}
        rlist = rlist.select(mask)
        assert len(rlist) > 0

        xc, yc, zc = rlist["xyzcal.px"].parts()
        xo, yo, zo = rlist["xyzobs.px.value"].parts()
        dx = xc - xo
        dy = yc - yo
        dz = zc - zo

        is_stills = dz.all_approx_equal(dz[0])

        histx = flex.histogram(dx, n_slots=100)
        histy = flex.histogram(dy, n_slots=100)
        Hxy, xedges, yedges = np.histogram2d(dx, dy, bins=(50, 50))

        if not is_stills:
            histz = flex.histogram(dz, n_slots=100)
            Hzy, zedges, yedges = np.histogram2d(dz, dy, bins=(50, 50))
            Hxz, xedges, zedges = np.histogram2d(dx, dz, bins=(50, 50))

        density_hist_layout = {
            "showlegend": False,
            "autosize": False,
            "width": 600,
            "height": 550,
            "margin": {"t": 50},
            "hovermode": "closest",
            "bargap": 0,
            "xaxis": {
                "domain": [0, 0.85],
                "showgrid": False,
                "zeroline": True,
                "zerolinewidth": 2,
                "zerolinecolor": "#969696",
                "title": "X",
            },
            "yaxis": {
                "domain": [0, 0.85],
                "showgrid": False,
                "zeroline": True,
                "zerolinewidth": 2,
                "zerolinecolor": "#969696",
                "title": "Y",
            },
            "xaxis2": {
                "domain": [0.85, 1],
                "showgrid": False,
                "zeroline": False,
                "bargap": 0,
            },
            "yaxis2": {
                "domain": [0.85, 1],
                "showgrid": False,
                "zeroline": False,
                "bargap": 0,
            },
        }

        d = {
            "residuals_xy": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": Hxy.transpose().tolist(),
                        "name": "density",
                        "colorscale": "Hot",
                        "reversescale": True,
                        "showscale": False,
                        "type": "contour",
                        "zsmooth": "best",
                    },
                    {
                        "x": list(histx.slot_centers()),
                        "y": list(histx.slots()),
                        "name": "dx histogram",
                        "marker": {"color": "rgb(102,0,0)"},
                        "yaxis": "y2",
                        "type": "bar",
                    },
                    {
                        "y": list(histy.slot_centers()),
                        "x": list(histy.slots()),
                        "name": "dy histogram",
                        "marker": {"color": "rgb(102,0,0)"},
                        "xaxis": "x2",
                        "type": "bar",
                        "orientation": "h",
                    },
                ],
                "layout": copy.deepcopy(density_hist_layout),
            }
        }
        d["residuals_xy"]["layout"]["title"] = "Centroid residuals in X and Y"

        if not is_stills:
            d["residuals_zy"] = {
                "data": [
                    {
                        "x": zedges.tolist(),
                        "y": yedges.tolist(),
                        "z": Hzy.transpose().tolist(),
                        "name": "density",
                        "colorscale": "Hot",
                        "reversescale": True,
                        "showscale": False,
                        "type": "contour",
                        "zsmooth": "best",
                    },
                    {
                        "x": list(histz.slot_centers()),
                        "y": list(histz.slots()),
                        "name": "dz histogram",
                        "marker": {"color": "rgb(102,0,0)"},
                        "yaxis": "y2",
                        "type": "bar",
                    },
                    {
                        "y": list(histy.slot_centers()),
                        "x": list(histy.slots()),
                        "name": "dy histogram",
                        "marker": {"color": "rgb(102,0,0)"},
                        "xaxis": "x2",
                        "type": "bar",
                        "orientation": "h",
                    },
                ],
                "layout": copy.deepcopy(density_hist_layout),
            }
            d["residuals_zy"]["layout"]["title"] = "Centroid residuals in Z and Y"
            d["residuals_zy"]["layout"]["xaxis"]["title"] = "Z"

            d["residuals_xz"] = {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": zedges.tolist(),
                        "z": Hxz.transpose().tolist(),
                        "name": "density",
                        "colorscale": "Hot",
                        "reversescale": True,
                        "showscale": False,
                        "type": "contour",
                        "zsmooth": "best",
                    },
                    {
                        "x": list(histx.slot_centers()),
                        "y": list(histx.slots()),
                        "name": "dx histogram",
                        "marker": {"color": "rgb(102,0,0)"},
                        "yaxis": "y2",
                        "type": "bar",
                    },
                    {
                        "y": list(histz.slot_centers()),
                        "x": list(histz.slots()),
                        "name": "dz histogram",
                        "marker": {"color": "rgb(102,0,0)"},
                        "xaxis": "x2",
                        "type": "bar",
                        "orientation": "h",
                    },
                ],
                "layout": copy.deepcopy(density_hist_layout),
            }
            d["residuals_xz"]["layout"]["title"] = "Centroid residuals in X and Z"
            d["residuals_xz"]["layout"]["yaxis"]["title"] = "Z"

        return d


class IntensityAnalyser:
    """Analyse the intensities."""

    def __init__(self, grid_size=None, pixels_per_bin=10):
        self.grid_size = grid_size
        self.pixels_per_bin = pixels_per_bin

        # Set the required fields
        self.required = ["intensity.sum.value", "intensity.sum.variance", "xyzcal.px"]

    def __call__(self, rlist):
        """Analyse the reflection centroids."""
        # FIXME Do the same and a comparison for intensity.prf

        # Check we have the required fields
        print("Analysing reflection intensities")
        if not ensure_required(rlist, self.required):
            return {"intensity": {}}

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
        Command.start(" Selecting only integrated reflections")
        mask = rlist.get_flags(rlist.flags.integrated, all=False)
        if mask.count(True) == 0:
            return {"intensity": {}}

        rlist = rlist.select(mask)
        Command.end(f" Selected {len(rlist)} integrated reflections")

        x, y, z = rlist["xyzcal.px"].parts()
        self.nbinsx, self.nbinsy = tuple(
            int(math.ceil(i))
            for i in (
                flex.max(x) / self.pixels_per_bin,
                flex.max(y) / self.pixels_per_bin,
            )
        )

        print(" Analysing distribution of I/Sigma")
        d = self.i_over_s_hist(rlist)
        print(" Analysing distribution of I/Sigma vs xy")
        d.update(self.i_over_s_vs_xy(rlist, "sum"))
        if "intensity.prf.value" in rlist:
            print(" Analysing distribution of I/Sigma vs xy")
            d.update(self.i_over_s_vs_xy(rlist, "prf"))
        print(" Analysing distribution of I/Sigma vs z")
        d.update(self.i_over_s_vs_z(rlist))
        print(" Analysing distribution of partialities")
        d.update(self.partiality_hist(rlist))
        print(" Analysing QE map")
        d.update(self.qe_xy(rlist))
        return {"intensity": d}

    def i_over_s_hist(self, rlist):
        """Analyse the correlations."""
        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        log_I_over_S = flex.log10(I_over_S)
        hist = flex.histogram(log_I_over_S, n_slots=20)

        return {
            "log_i_over_sigma_histogram": {
                "data": [
                    {
                        "x": list(hist.slot_centers()),
                        "y": list(hist.slots()),
                        "type": "bar",
                        "name": "log_i_over_sigma",
                    }
                ],
                "layout": {
                    "title": "Log I/σ(I) histogram",
                    "xaxis": {"title": "Log I/σ(I)"},
                    "yaxis": {"title": "Number of reflections"},
                    "bargap": 0,
                },
            }
        }

    def i_over_s_vs_xy(self, rlist, intensity_type):
        """Plot I/Sigma vs X/Y"""

        I_sig = flex.sqrt(rlist[f"intensity.{intensity_type}.variance"])
        I = rlist[f"intensity.{intensity_type}.value"]
        sel = (I_sig > 0) & (I > 0)
        rlist = rlist.select(sel)
        I = I.select(sel)
        I_sig = I_sig.select(sel)
        I_over_S = I / I_sig
        x, y, z = rlist["xyzcal.px"].parts()

        H, xedges, yedges = np.histogram2d(
            x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy)
        )
        H1, xedges, yedges = np.histogram2d(
            x.as_numpy_array(),
            y.as_numpy_array(),
            bins=(self.nbinsx, self.nbinsy),
            weights=flex.log10(I_over_S).as_numpy_array(),
        )

        nonzeros = np.nonzero(H)
        z = np.empty(H.shape)
        z[:] = np.NAN
        z[nonzeros] = H1[nonzeros] / H[nonzeros]

        return {
            f"i_over_sigma_{intensity_type}_vs_xy": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": z.transpose().tolist(),
                        "zmin": -1,
                        "zauto": False,
                        "type": "heatmap",
                        "name": f"i_over_sigma_{intensity_type}",
                        "colorbar": {"title": "Log I/σ(I)", "titleside": "right"},
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": f"Distribution of I({intensity_type})/σ vs X/Y",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            }
        }

    def i_over_s_vs_z(self, rlist):
        """Plot I/Sigma vs Z."""

        I = rlist["intensity.sum.value"]
        I_sig = flex.sqrt(rlist["intensity.sum.variance"])
        I_over_S = I / I_sig
        x, y, z = rlist["xyzcal.px"].parts()

        H, xedges, yedges = np.histogram2d(
            z.as_numpy_array(), flex.log10(I_over_S).as_numpy_array(), bins=(100, 100)
        )

        return {
            "i_over_sigma_vs_z": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": "i_over_sigma",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Distribution of I/σ(I) vs Z",
                    "xaxis": {"title": "Z", "showgrid": False},
                    "yaxis": {"title": "Log I/σ(I)", "showgrid": False},
                },
            }
        }

    def partiality_hist(self, rlist):
        """Analyse the partialities."""
        partiality = rlist["partiality"]
        hist = flex.histogram(partiality.select(partiality > 0), 0, 1, n_slots=20)

        return {
            "partiality_histogram": {
                "data": [
                    {
                        "x": list(hist.slot_centers()),
                        "y": list(hist.slots()),
                        "type": "bar",
                        "name": "partiality",
                    }
                ],
                "layout": {
                    "title": "Partiality histogram",
                    "xaxis": {"title": "Partiality", "range": [0, 1]},
                    "yaxis": {"title": "Number of reflections"},
                    "bargap": 0,
                },
            }
        }

    def qe_xy(self, rlist):
        """Look at the QE map in x, y"""

        if "qe" not in rlist:
            return {}

        assert len(rlist) > 0
        qe = rlist["qe"]

        xc, yc, zc = rlist["xyzcal.px"].parts()

        nbinsx, nbinsy = tuple(
            int(math.ceil(i))
            for i in (
                flex.max(xc) / self.pixels_per_bin,
                flex.max(yc) / self.pixels_per_bin,
            )
        )

        xc = xc.as_numpy_array()
        yc = yc.as_numpy_array()
        qe = qe.as_numpy_array()

        H, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy))
        H1, xedges, yedges = np.histogram2d(xc, yc, bins=(nbinsx, nbinsy), weights=qe)

        nonzeros = np.nonzero(H)
        z1 = np.empty(H.shape)
        z1[:] = np.NAN
        z1[nonzeros] = H1[nonzeros] / H[nonzeros]

        return {
            "qe_map": {
                "data": [
                    {
                        "name": "qe_map",
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": z1.transpose().tolist(),
                        "type": "heatmap",
                        "colorbar": {"title": "QE", "titleside": "right"},
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Calculated Quantum Efficiency",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            }
        }


class ZScoreAnalyser:
    """
    Analyse the distribution of intensity z-scores.

    z-scores are calculated as z = (I - <I>) / σ, where I is intensity,
    <I> is a weighted mean and σ is the measurement error in intensity,
    as calculated by dials.util.intensity_explorer.
    """

    def __init__(self):
        # Set the required fields
        self.required = [
            "miller_index",
            "intensity.sum.value",
            "intensity.sum.variance",
            "xyzobs.px.value",
        ]

    def __call__(self, rlist, experiments):
        """
        :param rlist: A reflection table containing, at the least, the following
        fields:
          * ``miller_index``
          * ``intensity.sum.value``
          * ``intensity.sum.variance``
          * ``xyzobs.px.value``
        :type rlist: `dials_array_family_flex_ext.reflection_table`
        :param experiments: An experiment list with space group information.
        :type experiments: `dxtbx.model.ExperimentList`
        :return: A dictionary describing several Plotly plots.
        :rtype:`dict`
        """

        from dials.util.intensity_explorer import IntensityDist

        # Check we have the required fields
        print("Analysing reflection intensities")
        if not ensure_required(rlist, self.required):
            return {}

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

        self.z_score_data = IntensityDist(rlist, experiments).rtable

        print(" Analysing distribution of intensity z-scores")
        d = self.z_score_hist()
        print(" Constructing normal probability plot of intensity z-scores")
        d.update(self.normal_probability_plot())
        print(" Plotting intensity z-scores against multiplicity")
        d.update(self.z_vs_multiplicity())
        print(" Plotting time series of intensity z-scores")
        d.update(self.z_time_series())
        print(" Plotting intensity z-scores versus weighted mean intensity")
        d.update(self.z_vs_I())
        print(" Plotting intensity z-scores versus I/sigma")
        d.update(self.z_vs_I_over_sigma())

        return d

    def z_score_hist(self):
        """
        Plot a histogram of z-scores.

        :return: A dictionary describing a Plotly plot.
        :rtype:`dict`
        """

        z_scores = self.z_score_data["intensity.z_score"]
        hist = flex.histogram(z_scores, -10, 10, 100)

        return {
            "z_score_histogram": {
                "data": [
                    {
                        "x": hist.slot_centers().as_numpy_array().tolist(),
                        "y": hist.slots().as_numpy_array().tolist(),
                        "type": "bar",
                        "name": "z_hist",
                    }
                ],
                "layout": {
                    "title": "Histogram of z-scores",
                    "xaxis": {"title": "z-score", "range": [-10, 10]},
                    "yaxis": {"title": "Number of reflections"},
                    "bargap": 0,
                },
            }
        }

    def normal_probability_plot(self):
        """
        Display a normal probability plot of the z-scores of the intensities.

        :return: A dictionary describing a Plotly plot.
        :rtype:`dict`
        """

        z_scores = self.z_score_data["intensity.z_score"]
        osm = self.z_score_data["intensity.order_statistic_medians"]

        return {
            "normal_probability_plot": {
                "data": [
                    {
                        "x": osm.as_numpy_array().tolist(),
                        "y": z_scores.as_numpy_array().tolist(),
                        "type": "scattergl",
                        "mode": "markers",
                        "name": "Data",
                    },
                    {
                        "x": [-5, 5],
                        "y": [-5, 5],
                        "type": "scatter",
                        "mode": "lines",
                        "name": "z = m",
                    },
                ],
                "layout": {
                    "title": "Normal probability plot",
                    "xaxis": {"title": "z-score order statistic medians, m"},
                    "yaxis": {
                        "title": "Ordered z-score responses, z",
                        "range": [-10, 10],
                    },
                },
            }
        }

    def z_vs_multiplicity(self):
        """
        Plot z-score as a function of multiplicity of observation.

        :return: A dictionary describing a Plotly plot.
        :rtype:`dict`
        """

        multiplicity = self.z_score_data["multiplicity"]
        z_scores = self.z_score_data["intensity.z_score"]

        return {
            "z_score_vs_multiplicity": {
                "data": [
                    {
                        "x": multiplicity.as_numpy_array().tolist(),
                        "y": z_scores.as_numpy_array().tolist(),
                        "type": "scattergl",
                        "mode": "markers",
                        "name": "Data",
                    }
                ],
                "layout": {
                    "title": "z-scores versus multiplicity",
                    "xaxis": {"title": "Multiplicity"},
                    "yaxis": {"title": "z-score"},
                },
            }
        }

    def z_time_series(self):
        """
        Plot a crude time-series of z-scores, with image number serving as a
        proxy for time.

        :return: A dictionary describing a Plotly plot.
        :rtype:`dict`
        """

        batch_number = self.z_score_data["xyzobs.px.value"].parts()[2]
        z_scores = self.z_score_data["intensity.z_score"]

        return {
            "z_score_time_series": {
                "data": [
                    {
                        "x": batch_number.as_numpy_array().tolist(),
                        "y": z_scores.as_numpy_array().tolist(),
                        "type": "scattergl",
                        "mode": "markers",
                        "name": "Data",
                    }
                ],
                "layout": {
                    "title": "z-scores versus image number",
                    "xaxis": {"title": "Image number"},
                    "yaxis": {"title": "z-score"},
                },
            }
        }

    def z_vs_I(self):
        """
        Plot z-scores versus intensity.

        :return: A dictionary describing a Plotly plot.
        :rtype:`dict`
        """

        intensity = self.z_score_data["intensity.mean.value"]
        z_scores = self.z_score_data["intensity.z_score"]

        return {
            "z_score_vs_I": {
                "data": [
                    {
                        "x": intensity.as_numpy_array().tolist(),
                        "y": z_scores.as_numpy_array().tolist(),
                        "type": "scattergl",
                        "mode": "markers",
                        "name": "Data",
                    }
                ],
                "layout": {
                    "title": "z-scores versus weighted mean intensity",
                    "xaxis": {"title": "Intensity (photon count)", "type": "log"},
                    "yaxis": {"title": "z-score"},
                },
            }
        }

    def z_vs_I_over_sigma(self):
        """
        Plot z-scores versus intensity.

        :return: A dictionary describing a Plotly plot.
        :rtype:`dict`
        """

        i_over_sigma = (
            self.z_score_data["intensity.mean.value"]
            / self.z_score_data["intensity.mean.std_error"]
        )
        z_scores = self.z_score_data["intensity.z_score"]

        return {
            "z_score_vs_I_over_sigma": {
                "data": [
                    {
                        "x": i_over_sigma.as_numpy_array().tolist(),
                        "y": z_scores.as_numpy_array().tolist(),
                        "type": "scattergl",
                        "mode": "markers",
                        "name": "Data",
                    }
                ],
                "layout": {
                    "title": "z-scores versus I/σ",
                    "xaxis": {"title": "I/σ", "type": "log"},
                    "yaxis": {"title": "z-score"},
                },
            }
        }


class ReferenceProfileAnalyser:
    """Analyse the reference profiles."""

    def __init__(self, grid_size=None, pixels_per_bin=10):
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
        """Analyse the reference profiles."""
        # Check we have the required fields
        print("Analysing reference profiles")
        if not ensure_required(rlist, self.required):
            return {"reference": {}}

        # Select only integrated reflections
        Command.start(" Selecting only integrated reflections")
        mask = rlist.get_flags(rlist.flags.integrated)
        if mask.count(True) == 0:
            return {"reference": {}}

        rlist = rlist.select(mask)
        Command.end(f" Selected {len(rlist)} integrated reflections")

        x, y, z = rlist["xyzcal.px"].parts()
        self.nbinsx, self.nbinsy = tuple(
            int(math.ceil(i))
            for i in (
                flex.max(x) / self.pixels_per_bin,
                flex.max(y) / self.pixels_per_bin,
            )
        )

        # Analyse distribution of reference spots
        print(" Analysing reference profile distribution vs x/y")
        d = self.reference_xy(rlist)
        print(" Analysing reference profile distribution vs z")
        d.update(self.reference_z(rlist))

        # Look at correlations between profiles
        def correlations(filename, rlist):
            """Call for reference spots and all reflections."""

            print(" Analysing reflection profile correlations")
            d = self.reflection_corr_hist(rlist, filename)

            print(" Analysing reflection profile correlations vs x/y")
            d.update(self.reflection_corr_vs_xy(rlist, filename))

            print(" Analysing reflection profile correlations vs z")
            d.update(self.reflection_corr_vs_z(rlist, filename))

            print(" Analysing reflection profile correlations vs I/Sigma")
            d.update(self.reflection_corr_vs_ios(rlist, filename))

            return d

        mask = rlist.get_flags(rlist.flags.reference_spot)
        d.update(self.reflection_correlations_vs_resolution(rlist))
        d.update(correlations("reference", rlist.select(mask)))
        d.update(correlations("reflection", rlist))

        return {"reference": d}

    def reflection_correlations_vs_resolution(self, rlist):
        """Analyse the distribution of reference profiles."""

        print(" Analysing reflection correlations vs resolution")
        from dials.algorithms.spot_finding.per_image_analysis import binner_d_star_cubed

        profile_correlation = rlist["profile.correlation"]
        d_spacings = rlist["d"]
        binner = binner_d_star_cubed(d_spacings)
        bin_centres = flex.double()
        ccs = flex.double()
        for d_bin in binner.bins:
            d_min = d_bin.d_min
            d_max = d_bin.d_max
            ds3_min = 1 / d_min**3
            ds3_max = 1 / d_max**3
            ds3_centre = (ds3_max - ds3_min) / 2 + ds3_min
            d_centre = 1 / ds3_centre ** (1 / 3)
            sel = (d_spacings < d_max) & (d_spacings >= d_min)
            if sel.count(True) == 0:
                continue
            bin_centres.append(d_centre)
            ccs.append(flex.mean(profile_correlation.select(sel)))

        d_star_sq_bins = uctbx.d_as_d_star_sq(bin_centres)

        def d_star_sq_to_d_ticks(d_star_sq, nticks):
            min_d_star_sq = min(d_star_sq)
            dstep = (max(d_star_sq) - min_d_star_sq) / nticks
            tickvals = [min_d_star_sq + (i * dstep) for i in range(nticks)]
            ticktext = [f"{uctbx.d_star_sq_as_d(dsq):.2f}" for dsq in tickvals]
            return tickvals, ticktext

        tickvals, ticktext = d_star_sq_to_d_ticks(d_star_sq_bins, nticks=5)

        return {
            "reflection_cc_vs_resolution": {
                "data": [
                    {
                        "x": list(d_star_sq_bins),  # d_star_sq
                        "y": list(ccs),
                        "type": "scatter",
                        "name": "profile_correlation",
                    }
                ],
                "layout": {
                    "title": "Reflection correlations vs resolution",
                    "xaxis": {
                        "title": "Resolution (Å)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {
                        "title": "Correlation with reference profile",
                        "range": [0, 1],
                    },
                },
            }
        }

    def reference_xy(self, rlist):
        """Analyse the distribution of reference profiles."""
        mask = rlist.get_flags(rlist.flags.reference_spot)
        rlist = rlist.select(mask)
        x, y, z = rlist["xyzcal.px"].parts()

        H, xedges, yedges = np.histogram2d(
            x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy)
        )

        return {
            "n_reference_profiles_vs_xy": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": "n_reference_profiles",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": "Reference profiles binned in X/Y",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            }
        }

    def reference_z(self, rlist):
        """Analyse the distribution of reference profiles."""
        x, y, z = rlist["xyzcal.px"].parts()
        hist = flex.histogram(z, n_slots=20)

        return {
            "n_reference_profiles_vs_z": {
                "data": [
                    {
                        "x": list(hist.slot_centers()),
                        "y": list(hist.slots()),
                        "type": "bar",
                        "name": "n_reference_profiles",
                    }
                ],
                "layout": {
                    "title": "Reference profiles binned in Z",
                    "xaxis": {"title": "Z"},
                    "yaxis": {"title": "Number of reflections"},
                    "bargap": 0,
                },
            }
        }

    def reflection_corr_hist(self, rlist, filename):
        """Analyse the correlations."""
        corr = rlist["profile.correlation"]
        hist = flex.histogram(corr, n_slots=20)

        return {
            f"{filename}_correlations_histogram": {
                "data": [
                    {
                        "x": list(hist.slot_centers()),
                        "y": list(hist.slots()),
                        "type": "bar",
                        "name": f"{filename}_correlations",
                    }
                ],
                "layout": {
                    "title": f"{filename.capitalize()} correlations histogram",
                    "xaxis": {"title": "Correlation with reference profile"},
                    "yaxis": {"title": "Number of reflections"},
                    "bargap": 0,
                },
            }
        }

    def reflection_corr_vs_xy(self, rlist, filename):
        """Analyse the correlations."""

        corr = rlist["profile.correlation"]
        x, y, z = rlist["xyzcal.px"].parts()

        H, xedges, yedges = np.histogram2d(
            x.as_numpy_array(), y.as_numpy_array(), bins=(self.nbinsx, self.nbinsy)
        )
        H1, xedges, yedges = np.histogram2d(
            x.as_numpy_array(),
            y.as_numpy_array(),
            bins=(self.nbinsx, self.nbinsy),
            weights=corr.as_numpy_array(),
        )

        nonzeros = np.nonzero(H)
        z = np.empty(H.shape)
        z[:] = np.NAN
        z[nonzeros] = H1[nonzeros] / H[nonzeros]

        return {
            f"{filename}_correlations_xy": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": z.transpose().tolist(),
                        "type": "heatmap",
                        "name": f"{filename}_correlations",
                        "colorbar": {
                            "title": "Correlation with reference profile",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                        "zmin": 0,
                        "zmax": 1,
                    }
                ],
                "layout": {
                    "title": f"{filename.capitalize()} correlations binned in X/Y",
                    "xaxis": {"domain": [0, 0.85], "title": "X", "showgrid": False},
                    "yaxis": {"title": "Y", "autorange": "reversed", "showgrid": False},
                    "width": 500,
                    "height": 450,
                },
            }
        }

    def reflection_corr_vs_z(self, rlist, filename):
        """Analyse the correlations."""

        corr = rlist["profile.correlation"]
        x, y, z = rlist["xyzcal.px"].parts()

        H, xedges, yedges = np.histogram2d(
            z.as_numpy_array(), corr.as_numpy_array(), bins=(100, 100)
        )

        return {
            f"{filename}_correlations_vs_z": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": f"{filename}_correlations",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": f"{filename.capitalize()} correlations vs Z",
                    "xaxis": {"title": "Z", "showgrid": False},
                    "yaxis": {
                        "title": "Correlation with reference profile",
                        "showgrid": False,
                    },
                },
            }
        }

    def reflection_corr_vs_ios(self, rlist, filename):
        """Analyse the correlations."""

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

        H, xedges, yedges = np.histogram2d(
            flex.log10(I_over_S).as_numpy_array(),
            corr.as_numpy_array(),
            bins=(100, 100),
        )

        return {
            f"{filename}_correlations_vs_ios": {
                "data": [
                    {
                        "x": xedges.tolist(),
                        "y": yedges.tolist(),
                        "z": H.transpose().tolist(),
                        "type": "heatmap",
                        "name": f"{filename}_correlations",
                        "colorbar": {
                            "title": "Number of reflections",
                            "titleside": "right",
                        },
                        "colorscale": "Viridis",
                    }
                ],
                "layout": {
                    "title": f"{filename.capitalize()} correlations vs Log I/σ(I)",
                    "xaxis": {"title": "Log I/σ(I)", "showgrid": False},
                    "yaxis": {
                        "title": "Correlation with reference profile",
                        "showgrid": False,
                    },
                },
            }
        }


def merging_stats_data(reflections, experiments):
    reflections["intensity"] = reflections["intensity.scale.value"]
    reflections["variance"] = reflections["intensity.scale.variance"]

    reflections = (
        reflections.split_by_experiment_id()
    )  # needed to get correct batch plots
    scaled_array = scaled_data_as_miller_array(reflections, experiments)
    stats_data = MergingStatisticsData(
        experiments, scaled_array, reflections=reflections
    )
    stats, anom_stats = merging_stats_from_scaled_array(
        scaled_array,
        n_bins=20,
        use_internal_variance=False,
    )
    stats_data.merging_statistics_result = stats
    stats_data.anom_merging_statistics_result = anom_stats
    stats_plots = make_merging_stats_plots(stats_data, run_xtriage_analysis=True)
    intensity_plots = stats_plots["anom_plots"]
    intensity_plots.update(
        i_over_sig_i_vs_i_plot(scaled_array.data(), scaled_array.sigmas())
    )
    return (
        stats_plots["resolution_plots"],
        stats_plots["xtriage_output"],
        stats_plots["misc_plots"],
        intensity_plots,
        stats_plots["batch_plots"],
        stats_plots["image_range_tables"],
        stats_plots["scaling_tables"],
    )


class Analyser:
    """Helper class to do all the analysis."""

    def __init__(self, params, grid_size=None, centroid_diff_max=1.5):
        """Setup the analysers."""
        self.params = params
        self.analysers = [
            StrongSpotsAnalyser(pixels_per_bin=self.params.pixels_per_bin),
            CentroidAnalyser(
                grid_size=grid_size,
                pixels_per_bin=self.params.pixels_per_bin,
                centroid_diff_max=centroid_diff_max,
            ),
            IntensityAnalyser(
                grid_size=grid_size, pixels_per_bin=self.params.pixels_per_bin
            ),
            ReferenceProfileAnalyser(
                grid_size=grid_size, pixels_per_bin=self.params.pixels_per_bin
            ),
        ]

    def __call__(self, rlist=None, experiments=None):
        """Do all the analysis."""
        json_data = {}

        if rlist is not None:
            for analyse in self.analysers:
                result = analyse(copy.deepcopy(rlist))
                if result is not None:
                    json_data.update(result)
        else:
            json_data.update(
                {"strong": {}, "centroid": {}, "intensity": {}, "reference": {}}
            )

        crystal_table = None
        expt_geom_table = None
        if experiments is not None:
            analyse = ScanVaryingCrystalAnalyser(self.params.orientation_decomposition)
            json_data.update(analyse(experiments))
            crystal_table, expt_geom_table = self.experiments_table(experiments)
            json_data.update(
                {
                    "scaling_model": {},
                    "xtriage_output": {},
                    "resolution_graphs": {},
                    "batch_graphs": {},
                    "misc_graphs": {},
                    "scaled_intensity_graphs": {},
                    "scaling_tables": [],
                    "image_range_tables": [],
                }
            )
            if any(experiments.scaling_models()):
                print("Analysing scaling model")
                d = {}
                for i, model in enumerate(experiments.scaling_models()):
                    if model is not None:
                        scaling_model_plots = plot_scaling_models(model)
                        if scaling_model_plots:
                            for name, plot in scaling_model_plots.items():
                                d.update({name + "_" + str(i): plot})
                json_data["scaling_model"] = d

            if rlist:
                if "inverse_scale_factor" in rlist:
                    print("Calculating and generating merging statistics plots")
                    try:
                        (
                            resolution_plots,
                            xtriage_output,
                            misc_plots,
                            scaled_intensity_plots,
                            batch_plots,
                            image_range_tables,
                            scaling_tables,
                        ) = merging_stats_data(rlist, experiments)
                    except DialsMergingStatisticsError as e:
                        print(f"Error merging stats data: {e}")
                    else:
                        json_data["resolution_graphs"] = resolution_plots
                        json_data["xtriage_output"] = xtriage_output
                        json_data["batch_graphs"] = batch_plots
                        json_data["misc_graphs"] = misc_plots
                        json_data["scaled_intensity_graphs"] = scaled_intensity_plots
                        json_data["scaling_tables"] = scaling_tables
                        json_data["image_range_tables"] = image_range_tables

        if self.params.output.html is not None:
            from jinja2 import ChoiceLoader, Environment, PackageLoader

            loader = ChoiceLoader(
                [
                    PackageLoader("dials", "templates"),
                    PackageLoader("dials", "static", encoding="utf-8"),
                ]
            )
            env = Environment(loader=loader)

            graphs = json_data

            import libtbx.load_env

            static_dir = libtbx.env.find_in_repositories("dials/static")

            if self.params.output.external_dependencies == "local":
                template = env.get_template("report_local_dep.html")
            elif self.params.output.external_dependencies == "embed":
                template = env.get_template("report_embed_dep.html")
            else:
                template = env.get_template("report.html")
            html = template.render(
                page_title="DIALS analysis report",
                scan_varying_graphs=graphs["scan_varying"],
                scaling_model_graphs=graphs["scaling_model"],
                xtriage_output=graphs["xtriage_output"],
                resolution_plots=graphs["resolution_graphs"],
                batch_graphs=graphs["batch_graphs"],
                image_range_tables=graphs["image_range_tables"],
                misc_plots=graphs["misc_graphs"],
                scaled_intensity_plots=graphs["scaled_intensity_graphs"],
                strong_graphs=graphs["strong"],
                centroid_graphs=graphs["centroid"],
                intensity_graphs=graphs["intensity"],
                reference_graphs=graphs["reference"],
                crystal_table=crystal_table,
                geometry_table=expt_geom_table,
                scaling_tables=graphs["scaling_tables"],
                static_dir=static_dir,
            )

            print(f"Writing html report to: {self.params.output.html}")
            with open(self.params.output.html, "wb") as f:
                f.write(html.encode("utf-8", "xmlcharrefreplace"))

        if self.params.output.json is not None:
            print(f"Writing json data to: {self.params.output.json}")
            with open(self.params.output.json, "w", encoding="utf-8") as f:
                json.dump(json_data, f)

    def experiments_table(self, experiments):
        assert experiments is not None

        crystal_table = []
        expt_geom_table = []

        latex_vector_template = (
            r"$$\left( \begin{array}{rrr} %5.4f & %5.4f & %5.4f \end{array} \right)$$"
        )
        latex_matrix_template = "\n".join(
            (
                r"$$\left( \begin{array}{rrr}",
                r"%5.4f & %5.4f & %5.4f \\",
                r"%5.4f & %5.4f & %5.4f \\",
                r"%5.4f & %5.4f & %5.4f \end{array} \right)$$",
            )
        )
        latex_unit_cell_template = r"$$\left( \begin{array}{cccccc} %5.3f & %5.3f & %5.3f & %5.3f & %5.3f & %5.3f \end{array}\right)$$"

        for expt in experiments:
            for panel_id, panel in enumerate(expt.detector):
                expt_geom_table.append(
                    (
                        "<strong>Panel %i</strong>:" % (panel_id + 1),
                        "Pixel size (mm):",
                        "{:.4f}, {:.4f}".format(*panel.get_pixel_size()),
                        "Image size (pixels):",
                        "%i, %i" % panel.get_image_size(),
                    )
                )
                expt_geom_table.append(
                    (
                        "",
                        "Trusted range:",
                        "{:g}, {:g}".format(*panel.get_trusted_range()),
                        "Thickness (mm):",
                        f"{panel.get_thickness():g}",
                    )
                )
                expt_geom_table.append(
                    (
                        "",
                        "Material:",
                        f"{panel.get_material()}",
                        "μ:",
                        f"{panel.get_mu():g}",
                    )
                )
                expt_geom_table.append(
                    (
                        "",
                        "Fast axis:",
                        latex_vector_template % panel.get_fast_axis(),
                        "Slow axis:",
                        latex_vector_template % panel.get_slow_axis(),
                    )
                )
                expt_geom_table.append(
                    (
                        "",
                        "Origin:",
                        latex_vector_template % panel.get_origin(),
                        "Distance (mm)",
                        f"{panel.get_distance():.4f}",
                    )
                )
                if len(expt.detector) == 1:
                    try:
                        # does the beam intersect with the panel?
                        panel.get_beam_centre(expt.beam.get_s0())
                    except RuntimeError:
                        continue
                    else:
                        expt_geom_table.append(
                            (
                                "",
                                "Max resolution (corners) (Å):",
                                "%.2f"
                                % panel.get_max_resolution_at_corners(
                                    expt.beam.get_s0()
                                ),
                                "Max resolution (inscribed circle) (Å):",
                                "%.2f"
                                % panel.get_max_resolution_ellipse(expt.beam.get_s0()),
                            )
                        )

            if expt.scan is not None:
                expt_geom_table.append(
                    (
                        "<strong>Scan:</strong>",
                        "Image range:",
                        "%i, %i" % expt.scan.get_image_range(),
                        "Oscillation:",
                        "{:.2f}&deg;, {:+.2f}&deg;/frame".format(
                            *expt.scan.get_oscillation()
                        ),
                    )
                )

            if expt.goniometer is not None:
                expt_geom_table.append(
                    (
                        "<strong>Goniometer:</strong>",
                        "Rotation axis:",
                        latex_vector_template % expt.goniometer.get_rotation_axis(),
                    )
                )
                expt_geom_table.append(
                    (
                        "",
                        "Fixed rotation:",
                        latex_matrix_template % expt.goniometer.get_fixed_rotation(),
                        "Setting rotation:",
                        latex_matrix_template % expt.goniometer.get_setting_rotation(),
                    )
                )

            if expt.crystal is not None:
                uc = expt.crystal.get_unit_cell().parameters()
                sgi = expt.crystal.get_space_group().info()
                umat = latex_matrix_template % expt.crystal.get_U()
                bmat = latex_matrix_template % expt.crystal.get_B()
                amat = latex_matrix_template % expt.crystal.get_A()
                crystal_table.append(
                    (
                        "<strong>Crystal:</strong>",
                        "Space group:",
                        sgi.symbol_and_number(),
                        "Unit cell:",
                        latex_unit_cell_template % uc,
                    )
                )
                crystal_table.append(
                    ("", "U matrix:", f"{umat}", "B matrix:", f"{bmat}")
                )
                crystal_table.append(("", "A = UB:", f"{amat}"))
                if expt.crystal.num_scan_points > 0:
                    abc = flex.vec3_double()
                    angles = flex.vec3_double()
                    for n in range(expt.crystal.num_scan_points):
                        (
                            a,
                            b,
                            c,
                            alpha,
                            beta,
                            gamma,
                        ) = expt.crystal.get_unit_cell_at_scan_point(n).parameters()
                        abc.append((a, b, c))
                        angles.append((alpha, beta, gamma))
                    a, b, c = abc.mean()
                    alpha, beta, gamma = angles.mean()
                    mean_uc = uctbx.unit_cell((a, b, c, alpha, beta, gamma))
                    crystal_table.append(
                        (
                            "",
                            "A sampled at %i scan points"
                            % expt.crystal.num_scan_points,
                            "",
                            "Average unit cell:",
                            latex_unit_cell_template % mean_uc.parameters(),
                        )
                    )

        # ensure all the rows are the same length
        for table in (expt_geom_table, crystal_table):
            for i_row in range(len(table)):
                while len(table[i_row]) < 5:
                    table[i_row] = list(table[i_row]) + [""]

        return crystal_table, expt_geom_table


class Script:
    """A class to encapsulate the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import ArgumentParser

        # Create the parser
        usage = "usage: dials.report [options] observations.refl"
        self.parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            read_reflections=True,
            read_experiments=True,
            check_format=False,
            epilog=help_message,
        )
        dials.util.log.print_banner()

    def run(self, args=None):
        """Run the script."""
        from dials.util.options import reflections_and_experiments_from_files

        # Parse the command line arguments
        params, options = self.parser.parse_args(args, show_diff_phil=True)

        # Show the help
        if len(params.input.reflections) != 1 and not len(params.input.experiments):
            self.parser.print_help()
            exit(0)

        reflections, experiments = reflections_and_experiments_from_files(
            params.input.reflections, params.input.experiments
        )

        # Analyse the reflections
        analyse = Analyser(
            params,
            grid_size=params.grid_size,
            centroid_diff_max=params.centroid_diff_max,
        )
        if len(reflections):
            reflections = reflections[0]
        else:
            reflections = None

        analyse(reflections, experiments)


@show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
