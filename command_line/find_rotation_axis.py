"""
Optimise the rotation axis orientation using the method from Gorelik et al.
(https://www.doi.org/10.1007/978-94-007-5580-2) with code adapted from
Stef Smeets's edtools (https://github.com/instamatic-dev/edtools).

Examples::

  dials.find_rotation_axis imported.expt strong.refl
"""

import logging
import sys

import libtbx.phil

import dials.util
import dials.util.log

# from dials.array_family import flex
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

# from dxtbx.model import ExperimentList


# Define a logger.
logger = logging.getLogger("dials.find_rotation_axis")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
xds_inp = None
    .type = path
    .help = "Path to XDS.INP file (also reads SPOT.XDS in the same directory)"

view = False
    .type = bool
    .help = "View phi/theta histogram with current rotation axis (omega)"

optimise = True
    .type = bool
    .help = "Optimise the rotation axis"

finetune = False
    .type = bool
    .help = "Fine-tune rotation axis from the value in XDS.INP or given with omega=VAL"

omega = None
    .type = float
    .help = "Use the given value of omega to plot the histogram or as starting point for the optimization"

opposite = False
    .type = bool
    .help = Try the opposite value as the one defined in XDS.INP (or as given by omega=VAL)"

output {
    experiments = optimised.expt
        .type = path
    log = "dials.find_rotation_axis.log"
        .type = str
}
"""
)

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def rotation_axis_to_xyz(rotation_axis, invert=False, setting="xds"):
    """Convert rotation axis angle to XYZ vector compatible with 'xds', or 'dials'
    Set invert to 'True' for anti-clockwise rotation
    """
    if invert:
        rotation_axis += np.pi

    rot_x = np.cos(rotation_axis)
    rot_y = np.sin(rotation_axis)
    rot_z = 0

    if setting == "dials":
        return rot_x, -rot_y, rot_z
    elif setting == "xds":
        return rot_x, rot_y, rot_z
    else:
        raise ValueError("Must be one of {'dials', 'xds'}")


def rotation_matrix(axis, theta):
    """Calculates the rotation matrix around axis of angle theta (radians)"""

    # axis = axis/np.sqrt(np.dot(axis,axis))

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
    """Take a set of reflections in XYZ and convert to polar (cylindrical) coordinates"""
    sx, sy, sz = arr.T
    out = np.empty((len(arr), 2))
    np.hypot(sx, sy, out=out[:, 0])
    np.arctan2(sz, out[:, 0], out=out[:, 1])
    np.arctan2(sy, sx, out=out[:, 0])
    return out


def cylinder_histo(xyz, bins=(1000, 500)):
    """Take reciprocal lattice vectors in XYZ format and output cylindrical projection.
    `Bins` gives the resolution of the 2D histogram."""
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


def plot_histo(H, xedges, yedges, title="Histogram"):
    """Plot the histogram of the cylindrical projection."""
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
    plt.show()


def make(arr, omega: float, wavelength: float):
    """
    Prepare xyz (reciprocal space coordinates) from reflection positions/angle (`arr`),
    which is the list of reflections read from XDS (SPOT.XDS)

    omega: rotation axis (degrees), which is defined by the angle between x
        (horizontal axis pointing right) and the rotation axis going in clockwise direction

    Note that:
        1. x<->y are flipped
    This is to ensure to match the XDS convention with the one I'm used to
    """

    reflections = arr[:, 0:2]
    angle = arr[:, 2]

    omega_rad = np.radians(omega)
    r = make_2d_rotmat(omega_rad)

    refs_ = np.dot(reflections, r)

    y, x = refs_.T  # NOTE 1

    R = 1 / wavelength
    C = R - np.sqrt(R ** 2 - x ** 2 - y ** 2).reshape(-1, 1)
    xyz = (
        np.c_[x * np.cos(angle), y, -x * np.sin(angle)]
        + C * np.c_[-np.sin(angle), np.zeros_like(angle), -np.cos(angle)]
    )

    return xyz


def optimize(
    arr,
    omega_start: float,
    wavelength=float,
    plusminus: int = 180,
    step: int = 10,
    hist_bins: (int, int) = (1000, 500),
    plot: bool = False,
) -> float:
    """
    Optimize the value of omega around the given point.

    omega_start: defines the starting angle
    step, plusminus: together with omega_start define the range of values to loop over
    hist_bins: size of the 2d histogram to produce the final phi/theta plot
    plot: toggle to plot the histogram after each step
    """

    r = np.arange(omega_start - plusminus, omega_start + plusminus, step)

    best_score = 0
    best_omega = 0

    for omega in r:
        xyz = make(arr, omega, wavelength)

        H, xedges, yedges = cylinder_histo(xyz, bins=hist_bins)

        var = np.var(H)

        logger.info(f"Omega: {omega:8.2f}, variance: {var:5.2f}")

        if plot:
            plot_histo(
                H,
                xedges,
                yedges,
                title=f"omega={omega:.2f}$^\\circ$ | variance: {var:.2f}",
            )

        xvals.append(omega)
        vvals.append(var)

        if var > best_score:
            best_omega = omega
            best_score = var

    logger.info(f"Best omega: {best_omega:.2f}; score: {best_score:.2f}")

    return best_omega


def parse_xds_inp(fn):
    """
    Parse the XDS.INP file to find the required numbers for the optimization
    Looks for wavelength, pixelsize, beam_center, oscillation range
    """
    with open(fn, "r") as f:
        for line in f:
            line = line.split("!", 1)[0].strip()
            match = False

            if "X-RAY_WAVELENGTH" in line:
                match = True
                wavelength = float(line.rsplit("X-RAY_WAVELENGTH=")[1].split()[0])
            if "ORGX=" in line:
                match = True
                orgx = float(line.rsplit("ORGX=")[1].split()[0])
            if "ORGY=" in line:
                match = True
                orgy = float(line.rsplit("ORGY=")[1].split()[0])
            if "OSCILLATION_RANGE=" in line:
                match = True
                osc_angle = float(line.rsplit("OSCILLATION_RANGE=")[1].split()[0])
            if "QX=" in line:
                match = True
                qx = float(line.rsplit("QX=")[1].split()[0])
            # if "QY=" in line:
            #    match = True
            #    qy = float(line.rsplit("QY=")[1].split()[0])
            if "DETECTOR_DISTANCE=" in line:
                match = True
                distance = float(line.rsplit("DETECTOR_DISTANCE=")[1].split()[0])
            if "ROTATION_AXIS=" in line:
                match = True
                inp = line.rsplit("ROTATION_AXIS=")[1].split()[0:3]
                rotx, roty, rotz = [float(val) for val in inp]

            if match:
                logger.info(line)

    omega_current = np.degrees(np.arctan2(roty, rotx))
    pixelsize = qx / (distance * wavelength)

    return np.array((orgx, orgy)), osc_angle, pixelsize, wavelength, omega_current


def load_spot_xds(fn, beam_center: [float, float], osc_angle: float, pixelsize: float):
    """
    Load the given SPOT.XDS file (`fn`) and return an array with the reciprocal
        x, y, and angle for the centroid of each reflection

    beam_center: coordinates of the primary beam, read from XDS.INP
    osc_angle: oscillation_angle (degrees) per frame, will be multiplied by the average frame number
        that a reflection appears on (column 3 in `arr`)
    pixelsize: defined in px/Ångström

    http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#SPOT.XDS
    """
    arr = np.loadtxt(fn)
    logger.info(arr.shape)

    osc_angle_rad = np.radians(osc_angle)

    reflections = arr[:, 0:2] - beam_center
    angle = arr[:, 2] * osc_angle_rad

    reflections *= pixelsize

    return np.c_[reflections, angle]


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

    parser = OptionParser(
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

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Configure the logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    # Check the models and data
    nexp = len(experiments)
    if nexp > 1:
        logger.info("Only the first experiment will be used for calculation")
    if nexp == 0 or len(reflections) == 0:
        parser.print_help()
        return
    if len(reflections) > 1:
        sys.exit("Only one reflections list can be imported at present")
    reflections = reflections[0]

    xds_inp = params.xds_inp
    if not xds_inp:
        xds_inp = Path("XDS.INP")
    else:
        xds_inp = Path(xds_inp)

    if not xds_inp.exists():
        logger.info(f"No such file: {xds_inp}\n")
        sys.exit()

    beam_center, osc_angle, pixelsize, wavelength, omega_current = parse_xds_inp(
        xds_inp
    )

    if params.omega is not None:
        omega_current = params.omega

    omega_opposite = omega_current + 180

    if params.opposite:
        omega_current = omega_opposite

    if omega_current > 180:
        omega_current -= 360

    if omega_opposite > 180:
        omega_opposite -= 360

    logger.info(f"\nBeam center: {beam_center[0]:.2f} {beam_center[1]:.2f}")
    logger.info(f"Oscillation angle (degrees): {osc_angle}")
    logger.info(f"Pixelsize: {pixelsize:.4f} px/Ångström")
    logger.info(f"Wavelength: {wavelength:.5f} Ångström")
    logger.info(f"Omega (current): {omega_current:.5f} degrees")
    logger.info(f"                 {np.radians(omega_current):.5f} radians")

    spot_xds = xds_inp.with_name("SPOT.XDS")

    if not spot_xds.exists():
        sys.exit(f"Cannot find file: {spot_xds}")

    arr = load_spot_xds(spot_xds, beam_center, osc_angle, pixelsize)

    hist_bins = 1000, 500

    if params.view:
        omega_final = omega_current
    elif params.optimise:
        global xvals
        global vvals
        xvals = []
        vvals = []

        omega_global = omega_local = omega_fine = 0

        if params.finetune:
            omega_tmp = omega_global = omega_current
        else:
            omega_tmp = 0
            omega_global = omega_tmp = optimize(
                arr, omega_tmp, wavelength, plusminus=180, step=5, hist_bins=hist_bins
            )

        omega_local = omega_tmp = optimize(
            arr, omega_tmp, wavelength, plusminus=5, step=1, hist_bins=hist_bins
        )

        omega_fine = omega_tmp = optimize(
            arr, omega_tmp, wavelength, plusminus=1, step=0.1, hist_bins=hist_bins
        )

        omega_final = omega_tmp

        logger.info("---")
        logger.info(f"Best omega (global search): {omega_global:.3f}")
        logger.info(f"Best omega (local search): {omega_local:.3f}")
        logger.info(f"Best omega (fine search): {omega_fine:.3f}")

    xyz = make(arr, omega_final, wavelength)
    H, xedges, yedges = cylinder_histo(xyz)

    var = np.var(H)
    logger.info(f"Variance: {var:.2f}")

    # check opposite
    xyz_opp = make(arr, omega_final + 180, wavelength)
    H_opp, xedges_opp, yedges_opp = cylinder_histo(xyz_opp)

    var_opp = np.var(H_opp)
    logger.info(f"Variance (opposite): {var_opp:.2f}")

    if var < var_opp:
        logger.info(
            f"\nOpposite angle ({omega_opposite:.2f} deg.) has higher variance!\n"
        )

    plot_histo(
        H, xedges, yedges, title=f"omega={omega_final:.2f}$^\\circ$ | var={var:.2f}"
    )

    if params.optimise and not params.view:
        # Plot rotation axis distribution curve
        plt.scatter(xvals, vvals, marker="+", lw=1.0, color="red")
        plt.xlabel("Rotation axis position ($^\\circ$)")
        plt.ylabel("Variance of the polar coordinate histogram")
        plt.title(f"Rotation axis determination | Maximum @ {omega_final:.2f}$^\\circ$")
        plt.show()

    omega_deg = omega_final
    omega_rad = np.radians(omega_final)

    logger.info(f"\nRotation axis found: {omega_deg:.2f} deg. / {omega_rad:.3f} rad.")

    logger.info(" - Instamatic (config/camera/camera_name.yaml)")
    omega_instamatic = omega_rad
    logger.info(f"    rotation_axis_vs_stage_xy: {omega_instamatic:.3f}")

    logger.info(" - XDS")
    rot_x_xds, rot_y_xds, rot_z_xds = rotation_axis_to_xyz(omega_rad, setting="xds")
    logger.info(f"    ROTATION_AXIS= {rot_x_xds:.4f} {rot_y_xds:.4f} {rot_z_xds:.4f}")
    logger.info(" - XDS (opposite rotation)")
    rot_x_xds, rot_y_xds, rot_z_xds = rotation_axis_to_xyz(
        omega_rad, setting="xds", invert=True
    )
    logger.info(f"    ROTATION_AXIS= {rot_x_xds:.4f} {rot_y_xds:.4f} {rot_z_xds:.4f}")

    logger.info(" - DIALS")
    rot_x_dials, rot_y_dials, rot_z_dials = rotation_axis_to_xyz(
        omega_rad, setting="dials"
    )
    logger.info(
        f"    geometry.goniometer.axes={rot_x_dials:.4f},{rot_y_dials:.4f},{rot_z_dials:.4f}"
    )
    logger.info(" - DIALS (opposite rotation)")
    rot_x_dials, rot_y_dials, rot_z_dials = rotation_axis_to_xyz(
        omega_rad, setting="dials", invert=True
    )
    logger.info(
        f"    geometry.goniometer.axes={rot_x_dials:.4f},{rot_y_dials:.4f},{rot_z_dials:.4f}"
    )

    logger.info(" - PETS (.pts)")
    omega_pets = omega_deg
    if omega_pets < 0:
        omega_pets += 360
    elif omega_pets > 360:
        omega_pets -= 360
    logger.info(f"    omega {omega_pets:.2f}")

    logger.info(" - RED (.ed3d)")
    omega_red = omega_deg
    if omega_red < -180:
        omega_red += 360
    elif omega_red > 180:
        omega_red -= 360
    logger.info(f"    ROTATIONAXIS    {omega_red:.4f}")


if __name__ == "__main__":
    run()
