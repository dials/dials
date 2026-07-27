from __future__ import annotations

from itertools import pairwise

from scitbx import matrix


class FrameOrientations:
    """Hold information about the experimental diffraction orientation
    (goniometer, beam, and crystal) at the centre of each frame of a
    rotation experiment."""

    def __init__(self, experiment):

        crystal = experiment.crystal
        beam = experiment.beam
        scan = experiment.scan
        gonio = experiment.goniometer

        image_range = scan.get_image_range()
        self.images = list(range(image_range[0], image_range[1] + 1))
        num_scan_points = scan.get_num_images() + 1

        # Extract experiment data at each scan point
        if beam.num_scan_points > 0:
            s0 = []
            for i in range(beam.num_scan_points):
                s0.append(matrix.col(beam.get_s0_at_scan_point(i)))
        else:
            s0 = [matrix.col(beam.get_s0()) for _ in range(num_scan_points)]

        if gonio.num_scan_points > 0:
            S_mats = [
                matrix.sqr(gonio.get_setting_rotation_at_scan_point(i))
                for i in range(gonio.num_scan_points)
            ]
        else:
            S_mats = [
                matrix.sqr(gonio.get_setting_rotation()) for _ in range(num_scan_points)
            ]

        F_mats = [
            matrix.sqr(gonio.get_fixed_rotation()) for _ in range(num_scan_points)
        ]
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

        # Sanity check that the number of scan points matches the number of images
        check = {len(x) for x in (s0, S_mats, F_mats, R_mats, U_mats)}
        assert len(check) == 1
        assert check.pop() == len(self.images) + 1

        # Store the beam vector at the frame centres by averaging between the scan points.
        self.s0 = []
        for d1, d2 in pairwise(s0):
            s0_frame = (d1 + d2) / 2
            self.s0.append(s0_frame)

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

        self.U = U_frames
        self.B = B_frames
        self.UB = [U * B for U, B in zip(U_frames, B_frames)]
