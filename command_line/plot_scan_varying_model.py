from __future__ import absolute_import, division, print_function

import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from dials.algorithms.refinement.rotation_decomposition import (
    solve_r3_rotation_for_angles_given_axes,
)
from dials.command_line.analyse_output import ensure_directory

import dials.util
from libtbx.phil import parse

phil_scope = parse(
    """
  output {
    directory = .
      .type = str
      .help = "The directory to store the results"

    format = *png pdf
      .type = choice

    debug = False
      .help = "print tables of values that will be plotted"
      .type = bool
      .expert_level = 1
  }

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

help_message = """

Generate plots of scan-varying models, including crystal orientation, unit cell
and beam centre, from the input refined.expt

Examples::

  dials.plot_scan_varying_model refined.expt

"""


class Script(object):
    """Class to run script."""

    def __init__(self):
        """Setup the script."""
        from dials.util.options import OptionParser

        usage = "usage: dials.plot_scan_varying_model [options] refined.expt"
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            check_format=False,
            epilog=help_message,
        )

    def run(self):
        """Run the script."""
        from dials.util.options import flatten_experiments
        from scitbx import matrix

        params, options = self.parser.parse_args()
        if len(params.input.experiments) == 0:
            self.parser.print_help()
            raise dials.util.Sorry("No experiments found in the input")
        experiments = flatten_experiments(params.input.experiments)

        # Determine output path
        self._directory = os.path.join(params.output.directory, "scan-varying_model")
        self._directory = os.path.abspath(self._directory)
        ensure_directory(self._directory)
        self._format = "." + params.output.format

        self._debug = params.output.debug

        # Decomposition axes
        self._e1 = params.orientation_decomposition.e1
        self._e2 = params.orientation_decomposition.e2
        self._e3 = params.orientation_decomposition.e3

        # cell plot
        dat = []
        for iexp, exp in enumerate(experiments):

            crystal = exp.crystal
            scan = exp.scan

            if crystal.num_scan_points == 0:
                print("Ignoring scan-static crystal")
                continue

            scan_pts = list(range(crystal.num_scan_points))
            cells = [crystal.get_unit_cell_at_scan_point(t) for t in scan_pts]
            cell_params = [e.parameters() for e in cells]
            a, b, c, aa, bb, cc = zip(*cell_params)
            phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
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
            try:
                cell_esds = [
                    crystal.get_cell_parameter_sd_at_scan_point(t) for t in scan_pts
                ]
                sig_a, sig_b, sig_c, sig_aa, sig_bb, sig_cc = zip(*cell_esds)
                cell_dat["sig_a"] = sig_a
                cell_dat["sig_b"] = sig_b
                cell_dat["sig_c"] = sig_c
                cell_dat["sig_aa"] = sig_aa
                cell_dat["sig_bb"] = sig_bb
                cell_dat["sig_cc"] = sig_cc
            except RuntimeError:
                pass

            if self._debug:
                print("Crystal in Experiment {}".format(iexp))
                print("Phi\ta\tb\tc\talpha\tbeta\tgamma\tVolume")
                msg = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}"
                line_dat = zip(phi, a, b, c, aa, bb, cc, vol)
                for line in line_dat:
                    print(msg.format(*line))
            dat.append(cell_dat)
        if dat:
            self.plot_cell(dat)

        # orientation plot
        dat = []
        for iexp, exp in enumerate(experiments):

            crystal = exp.crystal
            scan = exp.scan

            if crystal.num_scan_points == 0:
                print("Ignoring scan-static crystal")
                continue

            scan_pts = list(range(crystal.num_scan_points))
            phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
            Umats = [matrix.sqr(crystal.get_U_at_scan_point(t)) for t in scan_pts]
            if params.orientation_decomposition.relative_to_static_orientation:
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
                print("Crystal in Experiment {}".format(iexp))
                print("Image\tphi3\tphi2\tphi1")
                msg = "{0}\t{1}\t{2}\t{3}"
                line_dat = zip(phi, phi3, phi2, phi1)
                for line in line_dat:
                    print(msg.format(*line))
            dat.append(angle_dat)
        if dat:
            self.plot_orientation(dat)

        # beam centre plot
        dat = []
        for iexp, exp in enumerate(experiments):

            beam = exp.beam
            detector = exp.detector
            scan = exp.scan

            if beam.num_scan_points == 0:
                print("Ignoring scan-static beam")
                continue

            scan_pts = range(beam.num_scan_points)
            phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
            p = detector.get_panel_intersection(beam.get_s0())
            if p < 0:
                print("Beam does not intersect a panel")
                continue
            panel = detector[p]
            s0_scan_points = [
                beam.get_s0_at_scan_point(i) for i in range(beam.num_scan_points)
            ]
            bc_scan_points = [panel.get_beam_centre_px(s0) for s0 in s0_scan_points]
            bc_x, bc_y = zip(*bc_scan_points)
            dat.append({"phi": phi, "beam_centre_x": bc_x, "beam_centre_y": bc_y})
        if dat:
            self.plot_beam_centre(dat)

    def plot_cell(self, dat):
        plt.figure(figsize=(13, 10))
        gs = gridspec.GridSpec(4, 2, wspace=0.4, hspace=0.6)

        ax = plt.subplot(gs[0, 0])
        ax.ticklabel_format(useOffset=False)
        for cell in dat:
            if "sig_a" in cell:
                ax.errorbar(
                    cell["phi"][0::20], cell["a"][0::20], yerr=cell["sig_a"][0::20]
                )
            plt.plot(cell["phi"], cell["a"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"length $\left(\AA\right)$")
        plt.title("a")

        ax = plt.subplot(gs[0, 1])
        ax.ticklabel_format(useOffset=False)
        ymin, ymax = 0.0, 0.0
        for cell in dat:
            if "sig_aa" in cell:
                ax.errorbar(
                    cell["phi"][0::20], cell["alpha"][0::20], yerr=cell["sig_aa"][0::20]
                )
            plt.plot(cell["phi"], cell["alpha"])
            # choose the widest y range
            ymin = max(ymin, min(cell["alpha"]) - 0.05)
            ymax = max(ymax, max(cell["alpha"]) + 0.05)
            plt.axis(ymin=ymin, ymax=ymax)
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"$\alpha$")

        ax = plt.subplot(gs[1, 0])
        ax.ticklabel_format(useOffset=False)
        for cell in dat:
            if "sig_b" in cell:
                ax.errorbar(
                    cell["phi"][0::20], cell["b"][0::20], yerr=cell["sig_b"][0::20]
                )
            plt.plot(cell["phi"], cell["b"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"length $\left(\AA\right)$")
        plt.title("b")

        ax = plt.subplot(gs[1, 1])
        ax.ticklabel_format(useOffset=False)
        ymin, ymax = 0.0, 0.0
        for cell in dat:
            if "sig_bb" in cell:
                ax.errorbar(
                    cell["phi"][0::20], cell["beta"][0::20], yerr=cell["sig_bb"][0::20]
                )
            plt.plot(cell["phi"], cell["beta"])
            # choose the widest y range
            ymin = max(ymin, min(cell["beta"]) - 0.05)
            ymax = max(ymax, max(cell["beta"]) + 0.05)
            plt.axis(ymin=ymin, ymax=ymax)
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"$\beta$")

        ax = plt.subplot(gs[2, 0])
        ax.ticklabel_format(useOffset=False)
        for cell in dat:
            if "sig_c" in cell:
                ax.errorbar(
                    cell["phi"][0::20], cell["c"][0::20], yerr=cell["sig_c"][0::20]
                )
            plt.plot(cell["phi"], cell["c"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"length $\left(\AA\right)$")
        plt.title("c")

        ax = plt.subplot(gs[2, 1])
        ax.ticklabel_format(useOffset=False)
        ymin, ymax = 0.0, 0.0
        for cell in dat:
            if "sig_cc" in cell:
                ax.errorbar(
                    cell["phi"][0::20], cell["gamma"][0::20], yerr=cell["sig_cc"][0::20]
                )
            plt.plot(cell["phi"], cell["gamma"])
            # choose the widest y range
            ymin = max(ymin, min(cell["gamma"]) - 0.05)
            ymax = max(ymax, max(cell["gamma"]) + 0.05)
            plt.axis(ymin=ymin, ymax=ymax)
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"$\gamma$")

        ax = plt.subplot2grid((4, 2), (3, 0), colspan=2)
        ax.ticklabel_format(useOffset=False)
        for cell in dat:
            plt.plot(cell["phi"], cell["volume"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"volume $\left(\AA^3\right)$")
        plt.title("Cell volume")

        basename = os.path.join(self._directory, "unit_cell")
        fullname = basename + self._format
        print("Saving unit cell plot to {}".format(fullname))
        plt.savefig(fullname)

    def plot_orientation(self, dat):
        plt.figure(figsize=(13, 10))
        gs = gridspec.GridSpec(3, 1, wspace=0.4, hspace=0.6)

        ax = plt.subplot(gs[0, 0])
        ax.ticklabel_format(useOffset=False)
        for ori in dat:
            plt.plot(ori["phi"], ori["phi1"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"$\phi_1$")

        ax = plt.subplot(gs[1, 0])
        ax.ticklabel_format(useOffset=False)
        for ori in dat:
            plt.plot(ori["phi"], ori["phi2"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"$\phi_2$")

        ax = plt.subplot(gs[2, 0])
        ax.ticklabel_format(useOffset=False)
        for ori in dat:
            plt.plot(ori["phi"], ori["phi3"])
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"$\phi_3$")

        basename = os.path.join(self._directory, "orientation")
        fullname = basename + self._format
        print("Saving orientation plot to {}".format(fullname))
        plt.savefig(fullname)

    def plot_beam_centre(self, dat):
        plt.figure(figsize=(13, 10))
        gs = gridspec.GridSpec(2, 1, wspace=0.4, hspace=0.6)

        ax = plt.subplot(gs[0, 0])
        ax.ticklabel_format(useOffset=False)
        ymin, ymax = 0.0, 0.0
        for bc in dat:
            plt.plot(bc["phi"], bc["beam_centre_x"])
            ymin = max(ymin, min(bc["beam_centre_x"]) - 0.1)
            ymax = max(ymax, max(bc["beam_centre_x"]) + 0.1)
            plt.axis(ymin=ymin, ymax=ymax)
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"Beam centre X (pixels)")

        ax = plt.subplot(gs[1, 0])
        ax.ticklabel_format(useOffset=False)
        ymin, ymax = 0.0, 0.0
        for bc in dat:
            plt.plot(bc["phi"], bc["beam_centre_y"])
            ymin = max(ymin, min(bc["beam_centre_y"]) - 0.1)
            ymax = max(ymax, max(bc["beam_centre_y"]) + 0.1)
            plt.axis(ymin=ymin, ymax=ymax)
        plt.xlabel(r"rotation angle $\left(^\circ\right)$")
        plt.ylabel(r"angle $\left(^\circ\right)$")
        plt.title(r"Beam centre Y (pixels)")

        basename = os.path.join(self._directory, "beam_centre")
        fullname = basename + self._format
        print("Saving beam centre plot to {}".format(fullname))
        plt.savefig(fullname)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        script = Script()
        script.run()
