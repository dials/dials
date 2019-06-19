from __future__ import absolute_import, division, print_function

import os

from dxtbx.model import Crystal


def dump(experiments, directory):
    """
    Dump the experiments in mosflm format

    :param experiments: The experiments to dump
    :param directory: The directory to write to

    """
    from scitbx import matrix

    for i, experiment in enumerate(experiments):
        suffix = ""
        if len(experiments) > 1:
            suffix = "_%i" % (i + 1)

        sub_dir = "%s%s" % (directory, suffix)
        if not os.path.isdir(sub_dir):
            os.makedirs(sub_dir)
        detector = experiment.detector
        beam = experiment.beam
        scan = experiment.scan
        goniometer = experiment.goniometer

        # XXX imageset is getting the experimental geometry from the image files
        # rather than the input models.expt file
        imageset = experiment.imageset

        from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame

        R_to_mosflm = align_reference_frame(
            beam.get_s0(),
            (1.0, 0.0, 0.0),
            goniometer.get_rotation_axis(),
            (0.0, 0.0, 1.0),
        )

        cryst = experiment.crystal
        cryst = cryst.change_basis(
            cryst.get_space_group().info().change_of_basis_op_to_reference_setting()
        )
        A = matrix.sqr(cryst.get_A())
        A_inv = A.inverse()

        real_space_a = R_to_mosflm * A_inv.elems[:3]
        real_space_b = R_to_mosflm * A_inv.elems[3:6]
        real_space_c = R_to_mosflm * A_inv.elems[6:9]

        cryst_mosflm = Crystal(
            real_space_a,
            real_space_b,
            real_space_c,
            space_group=cryst.get_space_group(),
        )
        A_mosflm = matrix.sqr(cryst_mosflm.get_A())
        U_mosflm = matrix.sqr(cryst_mosflm.get_U())
        assert U_mosflm.is_r3_rotation_matrix(), U_mosflm
        w = beam.get_wavelength()

        index_mat = os.path.join(sub_dir, "index.mat")
        mosflm_in = os.path.join(sub_dir, "mosflm.in")
        print("Exporting experiment to %s and %s" % (index_mat, mosflm_in))

        with open(index_mat, "wb") as f:
            f.write(format_mosflm_mat(w * A_mosflm, U_mosflm, cryst.get_unit_cell()))

        img_dir, template = os.path.split(imageset.get_template())
        symmetry = cryst_mosflm.get_space_group().type().number()
        beam_centre = tuple(reversed(detector[0].get_beam_centre(beam.get_s0())))
        distance = detector[0].get_directed_distance()

        with open(mosflm_in, "wb") as f:
            f.write(
                write_mosflm_input(
                    directory=img_dir,
                    template=template,
                    symmetry=symmetry,
                    beam_centre=beam_centre,
                    distance=distance,
                    mat_file="index.mat",
                )
            )


def format_mosflm_mat(A, U, unit_cell, missets=(0, 0, 0)):
    lines = []
    uc_params = unit_cell.parameters()
    for i in range(3):
        lines.append(("%12.8f" * 3) % A.elems[i * 3 : 3 * (i + 1)])
    lines.append(("%12.3f" * 3) % missets)
    for i in range(3):
        lines.append("%12.8f" * 3 % U.elems[i * 3 : 3 * (i + 1)])
    lines.append(("%12.4f" * 6) % uc_params)
    lines.append(("%12.3f" * 3) % missets)
    return "\n".join(lines)


def write_mosflm_input(
    directory=None,
    template=None,
    symmetry=None,
    beam_centre=None,
    distance=None,
    mat_file=None,
):
    lines = []
    if directory is not None:
        lines.append("DIRECTORY %s" % directory)
    if template is not None:
        lines.append("TEMPLATE %s" % template)
    if symmetry is not None:
        lines.append("SYMMETRY %s" % symmetry)
    if beam_centre is not None:
        lines.append("BEAM %.3f %.3f" % beam_centre)
    if distance is not None:
        lines.append("DISTANCE %.4f" % distance)
    if mat_file is not None:
        lines.append("MATRIX %s" % mat_file)
    return "\n".join(lines)
