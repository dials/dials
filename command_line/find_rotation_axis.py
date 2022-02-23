"""
Optimise the rotation axis orientation using the method from Gorelik et al.
(https://doi.org/10.1007/978-94-007-5580-2_31) with code adapted from
Stef Smeets's edtools (https://github.com/instamatic-dev/edtools).

Three levels of search can be performed: a global search, a local search and a
fine search. At each trial position the sharpness of a cylindrical projection of
reciprocal lattice points is assessed via its variance. The correct orientation
maximises this variance.

Examples::

  dials.find_rotation_axis imported.expt strong.refl
"""

from __future__ import annotations

import logging
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import dxtbx.flumpy as flumpy
import libtbx.phil

import dials.util
import dials.util.log
from dials.array_family import flex

# from dials.array_family import flex
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

# Define a logger
logger = logging.getLogger("dials.find_rotation_axis")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
max_two_theta = 10.0
    .type = float
    .help = "Scattering angle limit to select reflections only in the central"
            "mostly flat region of the Ewald sphere surface"

view = False
    .type = bool
    .help = "View phi/theta histogram with current rotation axis (azimuth)"
            "instead of performing an optimisation"

global_search = True
    .type = bool
    .help = "Perform global search of the azimuthal angle. If False, only a"
            "local search will be performed around the start value or the value"
            "given by azimuth=VAL"

azimuth = None
    .type = float
    .help = "Use the given value of azimuth to plot the histogram or as the"
            "starting point for the optimization"

opposite = False
    .type = bool
    .help = "Try the opposite from the initial value (or that given by"
            "azimuth=VAL)"

output {
    experiments = optimised.expt
        .type = path

    log = "dials.find_rotation_axis.log"
        .type = str

    plot = find_rotation_axis
        .type = str
        .help = "Filename prefix for plots. Set to None for no plots"
}
"""
)


def rotation_axis_to_xyz(rotation_axis, invert=False):
    """
    Convert rotation axis angle to XYZ vector compatible with DIALS. Set invert
    to 'True' for anti-clockwise rotation
    """
    if invert:
        rotation_axis += np.pi

    rot_x = np.cos(rotation_axis)
    rot_y = np.sin(rotation_axis)
    rot_z = 0

    return rot_x, -rot_y, rot_z


def rotation_matrix(axis, theta):
    """Calculates the rotation matrix around axis of angle theta (radians)"""

    l = np.sqrt(np.dot(axis, axis))
    axis = axis / l

    a = np.cos(theta / 2)
    b, c, d = -1 * axis * np.sin(theta / 2)

    return np.array(
        [
            [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c],
        ]
    )


def make_2d_rotmat(theta):
    """Take angle in radians, and return 2D rotation matrix"""
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    return R


def random_sample(arr, n):
    """Select random sample of `n` rows from array"""
    indices = np.random.choice(arr.shape[0], n, replace=False)
    return arr[indices]


def xyz2cyl(arr):
    """
    Take a set of reflections in XYZ and convert to polar (cylindrical)
    coordinates
    """
    sx, sy, sz = arr.T
    out = np.empty((len(arr), 2))
    np.hypot(sx, sy, out=out[:, 0])
    np.arctan2(sz, out[:, 0], out=out[:, 1])
    np.arctan2(sy, sx, out=out[:, 0])
    return out


def cylinder_histo(xyz, bins=(1000, 500)):
    """
    Take reciprocal lattice vectors in XYZ format and output cylindrical
    projection. `bins` gives the resolution of the 2D histogram.
    """
    i, j = np.triu_indices(len(xyz), k=1)
    diffs = xyz[i] - xyz[j]
    polar = xyz2cyl(diffs)

    px, py = polar.T
    H, xedges, yedges = np.histogram2d(
        px,
        py,
        bins=bins,
        range=[
            [
                -np.pi,
                np.pi,
            ],
            [-np.pi / 2, np.pi / 2],
        ],
    )

    return H, xedges, yedges


def plot_histo(H, xedges, yedges, title="Histogram", filename=None):
    """Plot the histogram of the cylindrical projection."""
    plt.figure(figsize=(10, 7))
    plt.imshow(
        H.T,
        interpolation="nearest",
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        vmax=np.percentile(H, 99),
    )
    plt.title(title)
    plt.xlim(-np.pi, np.pi)
    plt.ylim(-np.pi / 2, np.pi / 2)
    plt.xlabel("phi ($\\pi$)")
    plt.ylabel("theta ($\\pi$)")
    if filename:
        plt.savefig(filename)
    else:
        plt.show()


def make(arr, azimuth: float, wavelength: float):
    """
    Prepare xyz (reciprocal space coordinates) from reflection positions/angle
    (`arr`). The azimuth angle (in degrees) gives the orientation of the
    rotation axis in the plane of the image, defined as the angle between x
    (horizontal axis pointing right) and the rotation axis going in a clockwise
    direction
    """

    reflections = arr[:, 0:2]
    angle = arr[:, 2]

    azimuth_rad = np.radians(azimuth)
    r = make_2d_rotmat(azimuth_rad)

    refs_ = np.dot(reflections, r)

    y, x = refs_.T

    R = 1 / wavelength
    C = R - np.sqrt(R**2 - x**2 - y**2).reshape(-1, 1)
    xyz = (
        np.c_[x * np.cos(angle), y, -x * np.sin(angle)]
        + C * np.c_[-np.sin(angle), np.zeros_like(angle), -np.cos(angle)]
    )

    return xyz


def optimise(
    arr,
    azimuth_start: float,
    wavelength=float,
    plusminus: int = 180,
    step: int = 10,
    hist_bins: (int, int) = (1000, 500),
    plot: bool = False,
) -> float:
    """
    Optimise the value of azimuth around the given point.

    Args:
        azimuth_start: defines the starting angle
        step, plusminus: together with azimuth_start define the range of values
            to loop over
        hist_bins: size of the 2d histogram to produce the final phi/theta plot
        plot: toggle to plot the histogram after each step

    Returns:
        The best value for the azimuth angle
    """

    r = np.arange(azimuth_start - plusminus, azimuth_start + plusminus, step)

    best_score = 0
    best_azimuth = 0

    for azimuth in r:
        xyz = make(arr, azimuth, wavelength)

        H, xedges, yedges = cylinder_histo(xyz, bins=hist_bins)

        var = np.var(H)

        logger.info(f"azimuth: {azimuth:8.2f}, variance: {var:5.2f}")

        if plot:
            plot_histo(
                H,
                xedges,
                yedges,
                title=f"azimuth={azimuth:.2f}$^\\circ$ | variance: {var:.2f}",
            )

        xvals.append(azimuth)
        vvals.append(var)

        if var > best_score:
            best_azimuth = azimuth
            best_score = var

    logger.info(f"Best azimuth: {best_azimuth:.2f}; score: {best_score:.2f}")

    return best_azimuth


def extract_spot_data(reflections, experiments, max_two_theta):
    """
    From the spot positions, extract reciprocal space X, Y and angle positions
    for each reflection up to the scattering angle max_two_theta
    """
    # Map reflections to reciprocal space
    reflections.centroid_px_to_mm(experiments)

    # Calculate scattering vectors
    reflections["s1"] = flex.vec3_double(len(reflections))
    reflections["2theta"] = flex.double(len(reflections))
    panel_numbers = flex.size_t(reflections["panel"])
    for i, expt in enumerate(experiments):
        if "imageset_id" in reflections:
            sel_expt = reflections["imageset_id"] == i
        else:
            sel_expt = reflections["id"] == i

        for i_panel in range(len(expt.detector)):
            sel = sel_expt & (panel_numbers == i_panel)
            x, y, _ = reflections["xyzobs.mm.value"].select(sel).parts()
            s1 = expt.detector[i_panel].get_lab_coord(flex.vec2_double(x, y))
            s1 = s1 / s1.norms() * (1 / expt.beam.get_wavelength())
            tt = s1.angle(expt.beam.get_s0(), deg=True)
            reflections["s1"].set_selected(sel, s1)
            reflections["2theta"].set_selected(sel, tt)

    # Filter reflections
    full_len = len(reflections)
    reflections = reflections.select(reflections["2theta"] <= max_two_theta)
    if len(reflections) < full_len:
        logger.info(
            f"{len(reflections)} reflections with 2θ ≤ {max_two_theta}° selected from {full_len} total"
        )

    x, y, _ = reflections["s1"].parts()
    _, _, angle = reflections["xyzobs.mm.value"].parts()
    arr = flumpy.to_numpy(x)
    arr = np.c_[
        x, flumpy.to_numpy(-y)
    ]  # Y is inverted to match calculation in edtools.find_rotation_axis
    arr = np.c_[arr, flumpy.to_numpy(angle)]
    return arr


@dials.util.show_mail_on_error()
def run(args=None, phil=phil_scope):
    """
    Load diffraction geometry and spot positions, optimise rotation axis
    orientation and save output files as specified.

    Args:
        args (list): Additional command-line arguments
        phil: The working PHIL parameters

    Returns:
        None
    """

    usage = "dials.find_rotation_axis [options] imported.expt strong.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the dials version
    logger.info(dials_version())

    # Log the difference between the PHIL scope definition and the active PHIL
    # scope, which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Check the models and data
    nexp = len(experiments)
    if nexp > 1:
        logger.info(
            "Only the first experiment will be used to determine oscillation and current rotation axis"
        )
        if len(experiments.goniometers()) > 1:
            sys.exit("Multiple experiments must share a goniometer model")
    if nexp == 0 or len(reflections) == 0:
        parser.print_help()
        return
    if len(reflections) > 1:
        sys.exit("Only one reflections list can be imported at present")
    reflections = reflections[0]
    expt = experiments[0]
    wavelength = expt.beam.get_wavelength()
    rotx, roty, _ = expt.goniometer.get_rotation_axis()
    azimuth_current = np.degrees(np.arctan2(roty, rotx))
    arr = extract_spot_data(reflections, experiments, params.max_two_theta)

    if params.azimuth is not None:
        azimuth_current = params.azimuth

    azimuth_opposite = azimuth_current + 180

    if params.opposite:
        azimuth_current = azimuth_opposite

    if azimuth_current > 180:
        azimuth_current -= 360

    if azimuth_opposite > 180:
        azimuth_opposite -= 360

    logger.info(
        f"Wavelength: {wavelength:.5f} Ångström\n"
        f"Azimuth (current): {azimuth_current:.5f} degrees\n"
        f"                 {np.radians(azimuth_current):.5f} radians"
    )

    hist_bins = 1000, 500

    if params.view:
        azimuth_final = azimuth_current
    else:
        global xvals
        global vvals
        xvals = []
        vvals = []

        azimuth_global = azimuth_local = azimuth_fine = 0

        if not params.global_search:
            azimuth_tmp = azimuth_global = azimuth_current
        else:
            logger.info("Performing global search")
            azimuth_tmp = 0
            azimuth_global = azimuth_tmp = optimise(
                arr, azimuth_tmp, wavelength, plusminus=180, step=5, hist_bins=hist_bins
            )

        logger.info("Performing local search")
        azimuth_local = azimuth_tmp = optimise(
            arr, azimuth_tmp, wavelength, plusminus=5, step=1, hist_bins=hist_bins
        )

        logger.info("Performing fine search")
        azimuth_fine = azimuth_tmp = optimise(
            arr, azimuth_tmp, wavelength, plusminus=1, step=0.1, hist_bins=hist_bins
        )

        azimuth_final = azimuth_tmp

        logger.info(
            "---\n"
            f"Best azimuth (global search): {azimuth_global:.3f}\n"
            f"Best azimuth (local search): {azimuth_local:.3f}\n"
            f"Best azimuth (fine search): {azimuth_fine:.3f}"
        )

    xyz = make(arr, azimuth_final, wavelength)
    H, xedges, yedges = cylinder_histo(xyz)

    var = np.var(H)
    logger.info(f"Variance: {var:.2f}")

    # check opposite
    xyz_opp = make(arr, azimuth_final + 180, wavelength)
    H_opp, xedges_opp, yedges_opp = cylinder_histo(xyz_opp)

    var_opp = np.var(H_opp)
    logger.info(f"Variance (opposite): {var_opp:.2f}")

    if var < var_opp:
        logger.info(
            f"\nOpposite angle ({azimuth_opposite:.2f}°) has higher variance!\n"
        )

    if params.output.plot:
        matplotlib.use("Agg")
        hist_filename = params.output.plot + "-histogram.png"
        plot_histo(
            H,
            xedges,
            yedges,
            title=f"azimuth={azimuth_final:.2f}$^\\circ$ | var={var:.2f}",
            filename=hist_filename,
        )

        if not params.view:
            # Plot rotation axis distribution curve
            plt.figure(figsize=(10, 7))
            plt.scatter(xvals, vvals, marker="+", lw=1.0, color="red")
            plt.xlabel("Rotation axis position ($^\\circ$)")
            plt.ylabel("Variance of the polar coordinate histogram")
            plt.title(
                f"Rotation axis determination | Maximum @ {azimuth_final:.2f}$^\\circ$"
            )
            projection_filename = params.output.plot + "-projection.png"
            plt.savefig(projection_filename)

    azimuth_deg = azimuth_final
    azimuth_rad = np.radians(azimuth_final)

    logger.info(
        f"\nRotation axis found: {azimuth_deg:.2f} deg. / {azimuth_rad:.3f} rad."
    )
    expt.goniometer.set_rotation_axis(rotation_axis_to_xyz(azimuth_rad))
    logger.info(str(expt.goniometer))

    logger.info(f"Saving optimised experiments to {params.output.experiments}")
    experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    run()
