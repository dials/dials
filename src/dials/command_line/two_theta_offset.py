from __future__ import annotations

from itertools import combinations

from libtbx.phil import parse
from scitbx import matrix
from scitbx.math import r3_rotation_axis_and_angle_from_matrix

import dials.util
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """

dials.two_theta_offset experiment_one.expt experiment_two.expt
"""

phil_scope = parse(
    """
offset_fast = 100.0
  .type = float
  .help = 'How far to move in the detector plane (fast direction)'
offset_slow = 100.0
  .type = float
  .help = 'How far to move in the detector plane (slow direction)'
min_distance = 10.0
  .type = float
  .help = 'Minimum shift in detector position'
""",
    process_includes=True,
)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        # The script usage
        usage = "usage: dials.two_theta_offset [options] experiment_one.expt experiment_two.expt"

        # Create the parser
        self.parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            check_format=False,
            read_experiments=True,
        )

    def run(self, args=None):
        """Execute the script."""
        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=True)

        # Check the number of experiments is at least 2
        experiments = flatten_experiments(params.input.experiments)
        if len(experiments) < 2:
            self.parser.print_help()
            return

        detectors = [experiment.detector[0] for experiment in experiments]

        for pair in combinations(detectors, 2):
            determine_axis(pair, params)

        crystals = [experiment.crystal for experiment in experiments]
        goniometers = [experiment.goniometer for experiment in experiments]

        FUs = []

        for c, g in zip(crystals, goniometers):
            u = matrix.sqr(c.get_U())
            f = matrix.sqr(g.get_fixed_rotation())
            FUs.append(f * u)

        for pair in combinations(FUs, 2):
            R = pair[1] * pair[0].inverse()
            rot = r3_rotation_axis_and_angle_from_matrix(R)
            angle = rot.angle(deg=True)
            axis = matrix.col(rot.axis)
            if abs(angle) < 10:
                continue
            print("Axis: %8.5f %8.5f %8.5f" % axis.elems, f"angle: {angle:7.4f}")


def determine_axis(detectors, params):
    offset_fast = params.offset_fast
    offset_slow = params.offset_slow
    min_distance = params.min_distance

    # pick two positions, at nominal origin offset in fast, slow

    x1 = matrix.col(detectors[0].get_origin())
    y1 = (
        matrix.col(detectors[0].get_origin())
        + offset_fast * matrix.col(detectors[0].get_fast_axis())
        + offset_slow * matrix.col(detectors[0].get_slow_axis())
    )

    x2 = matrix.col(detectors[1].get_origin())
    y2 = (
        matrix.col(detectors[1].get_origin())
        + offset_fast * matrix.col(detectors[1].get_fast_axis())
        + offset_slow * matrix.col(detectors[1].get_slow_axis())
    )

    # only allow this calculation if the detector has been moved a "significant"
    # amount

    if (x2 - x1).length() < min_distance:
        return

    centre, axis = find_centre_of_rotation(x1, x2, y1, y2)

    # compute "true" two-theta from these

    two_theta = component(x2 - centre, axis).angle(
        component(x1 - centre, axis), deg=True
    )

    print(
        "Centre: %7.4f %7.4f %7.4f" % centre.elems,
        "  axis: %7.4f %7.4f %7.4f" % axis.elems,
        f"angle: {two_theta:.2f}",
    )


def component(a, n):
    return a - a.dot(n) * n


def find_centre_of_rotation(x1, x2, y1, y2):
    """Find centre of rotation which takes position x1 -> x2 and y1 -> y2"""

    # chords of rotation of x, y

    cx = x2 - x1
    cy = y2 - y1

    # know axis is perpendicular to both of these -> is cross product

    axis = cx.cross(cy).normalize()

    # normal vector to y chord

    ny = component(cy, axis).normalize().cross(axis)

    # origin of normal vectors, centre of x, y chords

    ox = component(x1 + 0.5 * cx, axis)
    oy = component(y1 + 0.5 * cy, axis)

    # determine true origin of rotation - normal vector of x chord, construct
    # right-angle-triangle with hypotenuse from unknown origin of rotation
    # to central point of y chord oy, and adjacent the vector parallel to
    # reversed x chord => opposite is on vector from unknown origin of rotation
    # to ox

    ncx = cx.normalize()
    h = (oy - ox).dot(ncx)
    d = h / (ny).dot(-ncx)
    return oy + d * ny, axis


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
