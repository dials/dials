"""
Print a table of the orientation for every image of a dataset. The
orientation is expressed as a zone axis (a direction referenced to the direct
lattice) [uvw] giving the beam direction with respect to the crystal lattice.
Take into account any scan-varying models.

Usage: dials.frame_orientations refined.expt
"""


from __future__ import annotations

import sys
from itertools import tee

import matplotlib

from scitbx import matrix

import dials.util
from dials.util import tabulate
from dials.util.options import ArgumentParser, flatten_experiments

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# An itertools recipe in Python 3.7, but a module function in 3.10
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from libtbx.phil import parse

        # The phil scope
        phil_scope = parse(
            """
scale = unit *max_cell ewald_sphere_radius
    .type = choice
    .help = "Choose the scale for the direction vector in orthogonal"
            "coordinates prior to transformation into fractional"
            "coordinates [uvw]"

plot_filename = None
    .type = str
    .help = "Filename for a plot of angle between neighbouring frames"
            "(set to None for no plot)"
""",
            process_includes=True,
        )

        usage = "dials.frame_orientations refined.expt refined.refl"

        # Create the parser
        self.parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            check_format=False,
            epilog=__doc__,
        )

    def run(self, args=None):
        """Execute the script."""

        # Parse the command line
        self.params, _ = self.parser.parse_args(args, show_diff_phil=True)

        if not self.params.input.experiments:
            self.parser.print_help()
            sys.exit()

        # Try to load the models
        experiments = flatten_experiments(self.params.input.experiments)
        nexp = len(experiments)
        if nexp == 0:
            self.parser.print_help()
            sys.exit("No Experiments found in the input")

        # Set up a plot if requested
        if self.params.plot_filename:
            plt.figure()

        header = [
            "Image",
            "Beam direction (xyz)",
            "Zone axis [uvw]",
            "Angles between beam\nand axes a, b, c (deg)",
            "Angle from\nprevious image (deg)",
        ]
        for iexp, exp in enumerate(experiments):
            print(f"For Experiment id = {iexp}")
            print(exp.beam)
            print(exp.crystal)
            print(exp.scan)

            if self.params.scale == "ewald_sphere_radius":
                scale = 1.0 / exp.beam.get_wavelength()
            elif self.params.scale == "max_cell":
                uc = exp.crystal.get_unit_cell()
                scale = max(uc.parameters()[0:3])
            else:
                scale = 1.0
            print(
                "Beam direction scaled by {} = {:.3f} to "
                "calculate zone axis\n".format(self.params.scale, scale)
            )

            dat = extract_experiment_data(exp, scale)
            images = dat["images"]
            directions = dat["directions"]
            zone_axes = dat["zone_axes"]
            real_space_axes = dat["real_space_axes"]

            # calculate the angle between the beam and each crystal axis
            axis_angles = []
            for d, rsa in zip(directions, real_space_axes):
                angles = [d.angle(a, deg=True) for a in rsa]
                axis_angles.append("{:.2f} {:.2f} {:.2f}".format(*angles))

            # calculate the orientation offset between each image
            offset = [
                e1.angle(e2, deg=True) for e1, e2 in zip(zone_axes[:-1], zone_axes[1:])
            ]
            str_off = ["---"] + [f"{e:.8f}" for e in offset]

            rows = []
            for i, d, z, a, o in zip(
                images,
                directions,
                zone_axes,
                axis_angles,
                str_off,
            ):
                row = [
                    str(i),
                    "{:.8f} {:.8f} {:.8f}".format(*d.elems),
                    "{:.8f} {:.8f} {:.8f}".format(*z.elems),
                    a,
                    o,
                ]
                rows.append(row)

            # Print the table
            print(tabulate(rows, header))

            # Add to the plot, if requested
            if self.params.plot_filename:
                plt.scatter(images[1:], offset, s=1)

        # Finish and save plot, if requested
        if self.params.plot_filename:
            plt.xlabel("Image number")
            plt.ylabel(r"Angle from previous image $\left(^\circ\right)$")
            plt.title(r"Angle between neighbouring images")
            print(f"Saving plot to {self.params.plot_filename}")
            plt.savefig(self.params.plot_filename)

        print()


def extract_experiment_data(exp, scale=1):
    """Extract lists of the image number, beam direction and zone axis from an
    experiment"""
    crystal = exp.crystal
    beam = exp.beam
    scan = exp.scan
    gonio = exp.goniometer

    image_range = scan.get_image_range()
    images = list(range(image_range[0], image_range[1] + 1))
    num_scan_points = scan.get_num_images() + 1

    if beam.num_scan_points > 0:
        us0 = []
        for i in range(beam.num_scan_points):
            s0 = matrix.col(beam.get_s0_at_scan_point(i))
            us0.append(s0.normalize())
    else:
        us0 = [matrix.col(beam.get_unit_s0()) for _ in range(num_scan_points)]

    if gonio.num_scan_points > 0:
        S_mats = [
            matrix.sqr(gonio.get_setting_rotation_at_scan_point(i))
            for i in range(gonio.num_scan_points)
        ]
    else:
        S_mats = [
            matrix.sqr(gonio.get_setting_rotation()) for _ in range(num_scan_points)
        ]

    F_mats = [matrix.sqr(gonio.get_fixed_rotation()) for _ in range(num_scan_points)]
    start, stop = scan.get_array_range()
    R_mats = []
    axis = matrix.col(gonio.get_rotation_axis_datum())
    for i in range(start, stop + 1):
        phi = scan.get_angle_from_array_index(i, deg=False)
        R = matrix.sqr(axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=False))
        R_mats.append(R)

    if crystal.num_scan_points > 0:
        U_mats = [
            matrix.sqr(crystal.get_U_at_scan_point(i))
            for i in range(crystal.num_scan_points)
        ]
        B_mats = [
            matrix.sqr(crystal.get_B_at_scan_point(i))
            for i in range(crystal.num_scan_points)
        ]
    else:
        U_mats = [matrix.sqr(crystal.get_U()) for _ in range(num_scan_points)]
        B_mats = [matrix.sqr(crystal.get_B()) for _ in range(num_scan_points)]

    check = {len(x) for x in (us0, S_mats, F_mats, R_mats, U_mats)}
    assert len(check) == 1
    assert check.pop() == len(images) + 1

    # Construct full orientation matrix in the lab frame for each scan-point
    SRFU = (S * R * F * U for S, R, F, U in zip(S_mats, R_mats, F_mats, U_mats))

    # Now convert this to the orientation matrix at the centre of each frame by
    # calculating the linear transform that goes from the start of the frame
    # to the end, and then applying half of that to the start orientation
    U_frames = []
    for U1, U2 in pairwise(SRFU):
        M = U2 * U1.transpose()
        (
            angle,
            axis,
        ) = M.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
            deg=False
        )
        M_half = axis.axis_and_angle_as_r3_rotation_matrix(angle / 2, deg=False)
        U_frames.append(M_half * U1)

    # Convert the crystal B matrix at scan-points into a B matrix at the frame
    # centres. In this case, the transformation is not a rotation. Approximate
    # the answer by taking the average of the two B matrices.
    # FIXME is this kosher?
    B_frames = []
    for B1, B2 in pairwise(B_mats):
        B_frames.append((B1 + B2) / 2)

    UB_frames = [U * B for U, B in zip(U_frames, B_frames)]
    # This is the orthogonalisation matrix for the reciprocal space laboratory
    # frame. We want the real space fractionalisation matrix, which is its
    # transpose (https://dials.github.io/documentation/conventions.html)
    frac_mats = [m.transpose() for m in UB_frames]

    # Calculate zone axes, which also requires the beam directions at the frame
    # centres
    us0_frames = []
    for d1, d2 in pairwise(us0):
        us0_frames.append(((d1 + d2) / 2).normalize())
    zone_axes = [frac * (d * scale) for frac, d in zip(frac_mats, us0_frames)]

    # Now get the real space orthogonalisation matrix to calculate the real
    # space cell vectors at each image
    orthog_mats = (frac.inverse() for frac in frac_mats)
    h = matrix.col((1, 0, 0))
    k = matrix.col((0, 1, 0))
    l = matrix.col((0, 0, 1))
    real_space_axes = [(o * h, o * k, o * l) for o in orthog_mats]
    return {
        "images": images,
        "directions": us0_frames,
        "zone_axes": zone_axes,
        "real_space_axes": real_space_axes,
        "orientations": U_frames,
    }


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
