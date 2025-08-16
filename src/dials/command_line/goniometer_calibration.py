from __future__ import annotations

import iotbx.phil
from scitbx import matrix

import dials.util

help_message = """
dials.goniometer_calibration is a tool to aid calibration of multi-axis
goniometers.

The tool takes as input experiments.expt files for datasets recorded at the
goniometer datum setting and for each goniometer axis incremented in turn. It
outputs the axes and angles relating each consecutive pair of crystal setting
matrices in imgCIF and MOSFLM coordinate systems, and the CIF loop describing
the goniometer axes. Optionally it can also output an XOalign configuration file.

Either space_group must be specified or the parameter
use_space_group_from_experiments=True must be set

Examples::

dials.goniometer_calibration space_group=P422 \\
  experiments_o0_k0_p0.expt experiments_o0_k0_p48.expt \\
  experiments_o0_k48_p48.expt experiments_o48_k48_p48.expt
"""


phil_scope = iotbx.phil.parse(
    """
space_group = None
  .type = space_group
use_space_group_from_experiments = False
  .type = bool
output {
  xoalign = None
    .type = path
}
"""
)


@dials.util.show_mail_handle_errors()
def run(args=None):
    from dials.util.options import ArgumentParser, flatten_experiments

    usage = "dials.goniometer_calibration [options] models.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)
    if not params.use_space_group_from_experiments and params.space_group is None:
        parser.print_help()
        return

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) <= 1:
        parser.print_help()
        return

    from dials.algorithms.indexing.compare_orientation_matrices import (
        difference_rotation_matrix_axis_angle,
    )

    for experiment in experiments:
        crystal = experiment.crystal
        gonio = experiment.goniometer
        assert len(experiments) == (len(gonio.get_axes()) + 1)
        scan = experiment.scan
        fixed_rotation = matrix.sqr(gonio.get_fixed_rotation())
        setting_rotation = matrix.sqr(gonio.get_setting_rotation())
        rotation_axis = matrix.col(gonio.get_rotation_axis_datum())
        rotation_matrix = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
            scan.get_oscillation()[0], deg=True
        )
        U = matrix.sqr(crystal.get_U())
        U = setting_rotation * rotation_matrix * fixed_rotation * U
        crystal.set_U(U)
        if params.space_group is not None:
            crystal.set_space_group(params.space_group.group())

    rows = []

    from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame

    R_to_mosflm = align_reference_frame(
        experiments[0].beam.get_s0(),
        (1.0, 0.0, 0.0),
        experiments[0].goniometer.get_rotation_axis(),
        (0.0, 0.0, 1.0),
    )

    axes = []
    angles = []

    for i in range(len(experiments) - 1):
        target_angle = experiments[i + 1].goniometer.get_angles()[i]
        if i == experiments[i].goniometer.get_scan_axis():
            # rotation axis is canonical in our coordinate system
            axis = experiments[i].goniometer.get_axes()[i]
            angle = target_angle
        else:
            R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
                experiments[i].crystal,
                experiments[i + 1].crystal,
                target_angle=target_angle,
            )
        gonio = experiments[i + 1].goniometer
        axis_names = gonio.get_names()
        axes.append(axis)
        angles.append(angle)
        depends_on = "."
        if i + 1 < len(axis_names):
            depends_on = axis_names[i + 1]
        rows.insert(
            0,
            (
                axis_names[i],
                "rotation",
                "goniometer",
                depends_on,
                f"{axis[0]:.4f}",
                f"{axis[1]:.4f}",
                f"{axis[2]:.4f}",
                ".",
                ".",
                ".",
            ),
        )

    axis_names = experiments[0].goniometer.get_names()
    print("Goniometer axes and angles (ImgCIF coordinate system):")
    for axis, angle, name in zip(axes, angles, axis_names):
        print(
            f"{name}: ",
            f"rotation of {angle:.3f} degrees",
            "about axis ({:.5f},{:.5f},{:.5f})".format(*axis),
        )

    print()
    print("Goniometer axes and angles (MOSFLM coordinate system):")
    for axis, angle, name in zip(axes, angles, axis_names):
        print(
            f"{name}: ",
            f"rotation of {angle:.3f} degrees",
            "about axis ({:.5f},{:.5f},{:.5f})".format(
                *(R_to_mosflm * matrix.col(axis)).elems
            ),
        )

    print()
    print("ImgCIF _axis loop template:")
    from iotbx import cif

    loop = cif.model.loop(
        header=[
            "_axis.id",
            "_axis.type",
            "_axis.equipment",
            "_axis.depends_on",
            "_axis.vector[1]",
            "_axis.vector[2]",
            "_axis.vector[3]",
            "_axis.offset[1]",
            "_axis.offset[2]",
            "_axis.offset[3]",
        ]
    )
    for row in rows:
        loop.add_row(row)

    print(loop)

    if params.output.xoalign is not None:
        axes_mosflm = [(R_to_mosflm * matrix.col(axis)).elems for axis in axes]
        write_xoalign_config(
            params.output.xoalign, reversed(axes_mosflm), reversed(axis_names)
        )


def write_xoalign_config(file_name, axes, names):
    with open(file_name, "w") as f:
        print("GONIOMETER_AXES_NAMES = " + str(tuple(names)), file=f)
        print(
            "GONIOMETER_AXES = "
            + "["
            + ", ".join(
                "(" + ", ".join("%.5f" % x for x in axis) + ")" for axis in axes
            )
            + "]",
            file=f,
        )
        print("GONIOMETER_DATUM = (0,0,0) # in degrees", file=f)


if __name__ == "__main__":
    run()
