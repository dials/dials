from __future__ import annotations

import logging

import gemmi
import numpy as np

from dxtbx.model.experiment_list import ExperimentList
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from scitbx import matrix
from scitbx.math import r3_rotation_axis_and_angle_from_matrix

from dials.algorithms.refinement.rotation_decomposition import (
    solve_r3_rotation_for_angles_given_axes,
)
from dials.array_family import flex
from dials.command_line.frame_orientations import extract_experiment_data

logger = logging.getLogger(__name__)


def rotate_crystal(crystal, Rmat, axis, angle):

    Amats = []
    if crystal.num_scan_points > 0:
        scan_pts = list(range(crystal.num_scan_points))
        Amats = [
            Rmat
            * matrix.sqr(crystal.get_U_at_scan_point(t))
            * matrix.sqr(crystal.get_B_at_scan_point(t))
            for t in scan_pts
        ]

    crystal.rotate_around_origin(axis, angle, deg=False)
    if Amats:
        crystal.set_A_at_scan_points(Amats)

    return crystal


class PETSOutput:
    """Class to export integrated data in PETS CIF format"""

    def __init__(self, experiments, reflections, params):
        """Init with the parameters, experiments and reflections"""

        self.filename_prefix = params.pets.filename_prefix
        self.exp_id = params.pets.id
        self.partiality_cutoff = params.pets.partiality_cutoff
        self.flag_filtering = params.pets.flag_filtering
        self.excitation_error_cutoff = params.pets.virtual_frame.excitation_error_cutoff
        self.n_merged = params.pets.virtual_frame.n_merged
        self.step = params.pets.virtual_frame.step

        # Resolve a single intensity choice
        if "sum" in params.intensity and "profile" in params.intensity:
            logger.info(
                "For format=pets, intensity must be either 'sum' or 'profile'. Defaulting to 'sum'"
            )
            self.intensity_column = "intensity.sum.value"
            self.intensity_variance_column = "intensity.sum.variance"
        elif "sum" in params.intensity:
            if len(params.intensity) > 1:
                logger.info("Setting intensity=sum")
            self.intensity_column = "intensity.sum.value"
            self.intensity_variance_column = "intensity.sum.variance"
        elif "profile" in params.intensity:
            if len(params.intensity) > 1:
                logger.info("Setting intensity=profile")
            self.intensity_column = "intensity.prf.value"
            self.intensity_variance_column = "intensity.prf.variance"
        else:
            raise ValueError("The choice intensity= must be either 'sum' or 'profile'")

        self.experiments = experiments
        if len(reflections) != 1:
            raise ValueError("Only a single reflection file can be exported")
        self.reflections = reflections[0]

        self._check_experiments()
        self._check_reflections()
        self._reorient_coordinate_frame()

        # Calculate the frame orientation data
        self.frame_orientations = extract_experiment_data(self.experiment, scale=1)

        # Calculate the rotated missetting matrices. The orientations give the
        # full SRFMU matrix for each frame, where U is the scan-static
        # orientation matrix and M is the scan-varying misset. We want the
        # SRFM part, so post-multiply by the transpose of U
        Ut = matrix.sqr(self.experiment.crystal.get_U()).transpose()
        self.frame_orientations["rotated_missets"] = [
            M * Ut for M in self.frame_orientations["orientations"]
        ]

    @staticmethod
    def _rescale_zone_axis(axis):
        """Scale a vector so that the element with largest magnitude has
        magnitude 1.0, and also invert the vector. This is so the frame
        orientations zone axes match PETS format"""
        return -axis / max(axis.each_abs().elems)

    def _check_experiments(self):
        """Extract a single experiment from an experiment list and check that
        it is suitable for cif_pets output"""

        if self.exp_id is None:
            if len(self.experiments) != 1:
                raise ValueError("Please select a single experiment for export")
            else:
                self.exp_id = 0

        try:
            self.experiment = self.experiments[self.exp_id]
        except IndexError as e:
            raise IndexError(
                f"Unable to select experiment {self.exp_id} for export"
            ) from e

        # Check the unit cell is static
        crystal = self.experiment.crystal
        if crystal.num_scan_points > 0:
            Bmat = crystal.get_B()
            if not all(
                np.allclose(crystal.get_B_at_scan_point(i), Bmat)
                for i in range(crystal.num_scan_points)
            ):
                raise ValueError(
                    "The crystal must not have a scan-varying unit cell for export"
                )

        # Check no other models are scan-varying
        if (
            self.experiment.goniometer.num_scan_points > 0
            or self.experiment.beam.num_scan_points > 0
        ):
            raise ValueError("Only scan-varying crystal orientation is supported")

    def _check_reflections(self):
        """Check and filter the reflection table to include only the reflections
        for output"""

        self.reflections = self.reflections.select(
            self.reflections["id"] == self.exp_id
        )
        self.reflections["id"] *= 0

        required_keys = (self.intensity_column, self.intensity_variance_column)
        if any(x not in self.reflections for x in required_keys):
            raise ValueError(
                f"The reflection table requires these keys: {', '.join(required_keys)}"
            )

        # Select only fully-recorded reflections
        fulls = self.reflections["partiality"] >= self.partiality_cutoff
        logger.info(f"Removing {fulls.count(False)} partial reflections")
        self.reflections = self.reflections.select(fulls)

        # Select only "integrated" reflections, if requested
        if self.flag_filtering:
            if required_keys[0] == "intensity.sum.value":
                filter_flag = self.reflections.flags.integrated_sum
            else:
                filter_flag = self.reflections.flags.integrated_prf
            self.reflections = self.reflections.select(
                self.reflections.get_flags(filter_flag)
            )

    def _reorient_coordinate_frame(self):
        """Align a DIALS experiment and data in a reflection table to the
        PETS coordinate system. In that system, the s0 vector is aligned with
        -Z, while the rotation axis is in the X-Z plane, close to +X"""

        axis = matrix.col(self.experiment.goniometer.get_rotation_axis())
        us0 = matrix.col(self.experiment.beam.get_unit_s0())

        normal = us0.cross(-axis).normalize()
        orthogonalised_axis = us0.cross(normal).normalize()

        # Keep track of s1.us0 values to check reorientation
        s1_dot_us0 = self.reflections["s1"].dot(us0)

        R = align_reference_frame(us0, (0, 0, -1), orthogonalised_axis, (1, 0, 0))

        axis_angle = r3_rotation_axis_and_angle_from_matrix(R)
        axis = axis_angle.axis
        angle = axis_angle.angle(deg=False)
        logger.info(
            "Rotating experiment about axis ({:.4f}, {:.4f}, {:.4f}) by {:.3f}°".format(
                *axis, np.degrees(angle)
            )
        )

        self.experiment.detector.rotate_around_origin(axis, angle, deg=False)
        self.experiment.crystal = rotate_crystal(
            self.experiment.crystal, R, axis, angle
        )

        # Following does not work (https://github.com/cctbx/dxtbx/issues/454)
        # self.experiment.beam.rotate_around_origin(axis, angle, deg=False)
        # Set unit s0 and polarization normal directly instead (preserving inconsistent
        # beam, if that's what we have).
        new_us0 = (R * matrix.col(self.experiment.beam.get_unit_s0())).normalize()
        new_p_norm = (
            R * matrix.col(self.experiment.beam.get_polarization_normal())
        ).normalize()
        self.experiment.beam.set_unit_s0(new_us0)
        self.experiment.beam.set_polarization_normal(new_p_norm)

        # Rotating the goniometer is also complicated. See https://github.com/cctbx/dxtbx/pull/451
        new_datum = (
            R * matrix.col(self.experiment.goniometer.get_rotation_axis_datum())
        ).normalize()
        new_F = (
            R
            * matrix.sqr(self.experiment.goniometer.get_fixed_rotation())
            * R.transpose()
        )
        new_S = (
            R
            * matrix.sqr(self.experiment.goniometer.get_setting_rotation())
            * R.transpose()
        )
        self.experiment.goniometer.set_rotation_axis_datum(new_datum)
        self.experiment.goniometer.set_setting_rotation(new_S)
        self.experiment.goniometer.set_fixed_rotation(new_F)

        # Re-calculate s1 vectors and reciprocal lattice points with new geometry
        el = ExperimentList()
        el.append(self.experiment)
        self.reflections.map_centroids_to_reciprocal_space(el, calculated=True)

        error = flex.abs(s1_dot_us0 - self.reflections["s1"].dot(new_us0))
        if flex.max(error) > 1e-10:
            raise RuntimeError("Failed to rotate experiment correctly")

    def _set_virtual_frames(self):
        """Create a list of virtual frames for the experiment, each containing
        a table of reflections and the orientation information for the virtual
        frame"""

        # Get the frame orientation data
        # directions = self.frame_orientations["directions"]
        zone_axes = self.frame_orientations["zone_axes"]
        # real_space_axes = self.frame_orientations["real_space_axes"]
        rotated_missets = self.frame_orientations["rotated_missets"]

        # Get experiment geometry
        scan = self.experiment.scan
        s0 = matrix.col(self.experiment.beam.get_s0())
        inv_wl = s0.length()
        A = matrix.sqr(self.experiment.crystal.get_A())

        _, _, frames = self.reflections["xyzcal.px"].parts()
        scan = self.experiment.scan
        self.virtual_frames = []
        arr_start, arr_end = scan.get_array_range()
        # Loop over the virtual frames
        starts = list(range(arr_start, arr_end, self.step))
        for start in starts:
            stop = start + self.n_merged
            if stop > arr_end:
                break

            # Look up the orientation data using an index, which is the centre
            # of the virtual frame, offset so that the scan starts from 0
            centre = (start + stop) / 2.0
            index = int(centre) - arr_start
            M = rotated_missets[index]
            za = zone_axes[index]

            # If the number of merged frames is even, then we need to average
            # this orientation data with that from the previous frame to get the
            # orientation in the centre of the virtual frame
            if centre - int(centre) == 0.0:
                Mprev = rotated_missets[index - 1]
                RR = M * Mprev.transpose()
                (
                    angle,
                    axis,
                ) = RR.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
                    deg=False
                )
                R = axis.axis_and_angle_as_r3_rotation_matrix(angle / 2, deg=False)
                M = R * Mprev

                za_prev = zone_axes[index - 1]
                delta = za - za_prev
                za = za_prev + delta / 2

            # Rescale zone axis to match PETS format
            za = self._rescale_zone_axis(za)

            # Calculate the excitation error at this orientation
            UB = flex.mat3_double(len(self.reflections), M * A)
            r = UB * self.reflections["miller_index"].as_vec3_double()
            s1 = r + s0
            excitation_err = inv_wl - s1.norms()

            # Select only the reflections within the excitation error cutoff
            # and sort
            refs = self.reflections.select(
                abs(excitation_err) <= self.excitation_error_cutoff
            )
            refs.sort("miller_index")

            self.virtual_frames.append(
                {
                    "reflections": refs,
                    "rotated_misset": M,
                    "zone_axis": za,
                }
            )

    def write_dyn_cif_pets(self):

        self._set_virtual_frames()

        cif_filename = self.filename_prefix + ".cif_pets"
        cell_params = self.experiment.crystal.get_unit_cell().parameters()
        volume = self.experiment.crystal.get_unit_cell().volume()
        wavelength = self.experiment.beam.get_wavelength()
        UBmatrix = matrix.sqr(self.experiment.crystal.get_A()).as_numpy_array()
        # tilt semi-angle of the virtual frame
        precession_angle = (
            self.experiment.scan.get_oscillation()[1] * self.n_merged
        ) / 2

        # Some of these values are left as dummy zeroes for DIALS
        comment = f""";
data collection geometry: continuous rotation
dstarmax:  0.000
RC width:  0.00000
mosaicity:  0.000
rotation axis position:  000.000
reflection size:    0.000
Virtual frame settings: number of merged frames:  {self.n_merged}
                        step between frames:      {self.step}
                        sum all intensities:      1
;"""
        uvws = []
        for frame_id, virtual_frame in enumerate(self.virtual_frames):
            u, v, w = virtual_frame["zone_axis"]
            R = virtual_frame["rotated_misset"]
            # Decompose U = Rω * Rα * Rβ, where:
            # α is around 1,0,0
            # β is around 0,1,0
            # ω is around 0,0,1
            omega, alpha, beta = solve_r3_rotation_for_angles_given_axes(
                R, (0, 0, 1), (1, 0, 0), (0, 1, 0), deg=True
            )
            scale = 1  # dummy value
            uvws += [
                [frame_id + 1, u, v, w, precession_angle, alpha, beta, omega, scale]
            ]

        ref_list = []
        for frame_id, virtual_frame in enumerate(self.virtual_frames):
            refs = virtual_frame["reflections"]
            refs["intensity.sigma"] = flex.sqrt(refs[self.intensity_variance_column])
            for r in refs.rows():
                h, k, l = r["miller_index"]
                i = r[self.intensity_column]
                sig_i = r["intensity.sigma"]
                ref_list.append([h, k, l, i, sig_i, frame_id + 1])

        self._write_dyn_cif_pets(
            cif_filename,
            cell_params,
            volume,
            wavelength,
            UBmatrix,
            uvws,
            ref_list,
            comment,
        )

    def _write_dyn_cif_pets(
        self, name, lat_params, vol, lam, UBmatrix, uvws, ref_list, comment
    ):
        """
        - name : name of dataset
        - lat_params : (a,b,c,alpha,beta,gamma)
        - vol : volume
        - lam : wavelength
        - UBmatrix : 3x3 matrix orientation
        - uvws : list of lists with elements corresponding to the diffrn_zone_axis loop
        - ref_list : list of lists with elements corresponding to the refln_ loop
        """
        doc = gemmi.cif.Document()
        b = doc.add_new_block("pets")

        pair = [
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            "_cell_volume",
            "_diffrn_radiation_wavelength",
            "_diffrn_orient_matrix_UB_11",
            "_diffrn_orient_matrix_UB_12",
            "_diffrn_orient_matrix_UB_13",
            "_diffrn_orient_matrix_UB_21",
            "_diffrn_orient_matrix_UB_22",
            "_diffrn_orient_matrix_UB_23",
            "_diffrn_orient_matrix_UB_31",
            "_diffrn_orient_matrix_UB_32",
            "_diffrn_orient_matrix_UB_33",
            "_diffrn_pets_omega",
        ]
        values = np.hstack([lat_params, [vol, lam], UBmatrix.flatten(), [0]])
        for k, v in zip(pair, values):
            b.set_pair(k, "%.5f" % v)
        b.set_pair("_diffrn_measurement_details", comment)

        #### orientation info
        cols = [
            "id",
            "u",
            "v",
            "w",
            "precession_angle",
            "alpha",
            "beta",
            "omega",
            "scale",
        ]
        lo = b.init_loop("_diffrn_zone_axis_", cols)
        row_fmt = "{:4}," + "{:8.5f}," * 3 + "{:12.3f}," + "{:9.3f}," * 3 + "{:9.3f}"
        for uvw in uvws:
            lo.add_row(row_fmt.format(*uvw).split(","))

        #### Intensities
        logger.info("Writing intensities to cif")
        colsI = [
            "index_h",
            "index_k",
            "index_l",
            "intensity_meas",
            "intensity_sigma",
            "zone_axis_id",
        ]
        li = b.init_loop("_refln_", colsI)
        row_fmt = " " + "{:3}," * 3 + "{:11.2f}," * 2 + "{:5}"
        for row in ref_list:
            li.add_row(row_fmt.format(*row).split(","))

        doc.write_file(name)
        logger.info(f"File {name} written")
        return
