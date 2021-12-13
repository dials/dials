import logging

import gemmi
import numpy as np
import pandas as pd

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

        self.filename_prefix = params.filename_prefix
        self.exp_id = params.id
        self.partiality_cutoff = params.partiality_cutoff
        self.excitation_error_cutoff = params.virtual_frame.excitation_error_cutoff
        self.n_merged = params.virtual_frame.n_merged
        self.step = params.virtual_frame.step
        self.experiments = experiments
        if len(reflections) != 1:
            raise ValueError("Only a single reflection file can be exported")
        self.reflections = reflections[0]

        self._check_experiments()
        self._check_reflections()
        self._reorient_coordinate_frame()

        # Calculate the frame orientation data
        self.frame_orientations = extract_experiment_data(self.experiment, scale=1)

        # Calculate the frame orientation data (in the DIALS coordinate frame)
        self.frame_orientations = extract_experiment_data(self.experiment, scale=1)

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

        required_keys = ("intensity.sum.value", "intensity.sum.variance")
        if any(x not in self.reflections for x in required_keys):
            raise ValueError(
                f"The reflection table requires these keys: {', '.join(required_keys)}"
            )

        # Select only fully-recorded reflections
        fulls = self.reflections["partiality"] >= self.partiality_cutoff
        logger.info(f"Removing {fulls.count(False)} partial reflections")
        self.reflections = self.reflections.select(fulls)

        # Select only integrated reflections
        self.reflections = self.reflections.select(
            self.reflections.get_flags(self.reflections.flags.integrated)
        )

    def _reorient_coordinate_frame(self):
        """Align a DIALS experiment and data in a reflection table to the
        PETS coordinate system. In that system, the s0 vector is aligned with
        -Z, while the rotation axis is in the X-Z plane, close to +X"""

        axis = matrix.col(self.experiment.goniometer.get_rotation_axis())
        us0 = matrix.col(self.experiment.beam.get_unit_s0())

        normal = us0.cross(-axis).normalize()
        orthogonalised_axis = us0.cross(normal).normalize()

        R = align_reference_frame(us0, (0, 0, -1), orthogonalised_axis, (1, 0, 0))

        axis_angle = r3_rotation_axis_and_angle_from_matrix(R)
        axis = axis_angle.axis
        angle = axis_angle.angle()
        logger.info(f"Rotating experiment about axis {axis} by {np.degrees(angle)}°")

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

    def _set_virtual_frames(self):
        """Create a list of virtual frames for the experiment, each containing
        a table of reflections and the orientation information for the virtual
        frame"""

        # Get the frame orientation data
        # directions = self.frame_orientations["directions"]
        zone_axes = self.frame_orientations["zone_axes"]
        # real_space_axes = self.frame_orientations["real_space_axes"]
        orientations = self.frame_orientations["orientations"]

        # Get experiment geometry
        scan = self.experiment.scan
        axis = matrix.col(self.experiment.goniometer.get_rotation_axis())
        s0 = matrix.col(self.experiment.beam.get_s0())
        inv_wl = s0.length()

        _, _, frames = self.reflections["xyzcal.px"].parts()
        scan = self.experiment.scan
        self.virtual_frames = []
        arr_start, arr_end = scan.get_array_range()
        # Loop over the virtual frames.
        starts = list(range(arr_start, arr_end, self.step))
        for start in starts:
            stop = start + self.n_merged
            if stop > arr_end:
                stop = arr_end

            # Take only reflections whose centroids are inside this virtual frame
            sel = (frames > start) & (frames <= stop)
            refs = self.reflections.select(sel)

            centre = (start + stop) / 2.0
            # alpha_start = scan.get_angle_from_array_index(start, deg=False)
            # alpha_stop = scan.get_angle_from_array_index(stop, deg=False)
            alpha_centre = scan.get_angle_from_array_index(centre, deg=False)

            # The relps are given at zero rotation angle. In order to calculate
            # extinction distance, it's quickest to rotate the Ewald sphere
            s0_centre = s0.rotate_around_origin(axis, -alpha_centre, deg=False)
            rotated_s1 = refs["rlp"] + s0_centre
            extinction = rotated_s1.norms() - inv_wl

            # Now select only the reflections within the extinction distance
            # cutoff.
            refs.select(abs(extinction) <= self.excitation_error_cutoff)

            # Look up the orientation data using an index, which is the centre
            # of the virtual frame, offset so that the scan starts from 0
            index = int(centre) - arr_start

            self.virtual_frames.append(
                {
                    "reflections": refs,
                    "orientation": orientations[index],
                    "zone_axis": zone_axes[index],
                }
            )

    def write_cif_pets(self):
        pass

    def write_dyn_cif_pets(self):

        self._set_virtual_frames()

        cif_filename = self.filename_prefix + "_dyn.cif_pets"
        (
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        ) = self.experiment.crystal.get_unit_cell().parameters()
        volume = self.experiment.crystal.get_unit_cell().volume()
        wavelength = self.experiment.beam.get_wavelength()
        UBmatrix = matrix.sqr(self.experiment.crystal.get_A()).as_numpy_array()

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
            precession_angle = 0.0  # dummy value
            U = virtual_frame["orientation"]
            # Decompose U = Rω * Rβ * Rα, where:
            # α is around 1,0,0
            # β is around 0,1,0
            # ω is around 0,0,1
            # FIXME this is not confirmed as matching PETS2 yet
            alpha, beta, omega = solve_r3_rotation_for_angles_given_axes(
                U, (0, 0, 1), (0, 1, 0), (1, 0, 0), deg=True
            )
            scale = 1  # dummy value
            uvws += [
                [frame_id + 1, u, v, w, precession_angle, alpha, beta, omega, scale]
            ]

        dfI = pd.DataFrame(
            columns=[
                "index_h",
                "index_k",
                "index_l",
                "intensity_meas",
                "intensity_sigma",
                "F",
            ]
        )
        for frame_id, virtual_frame in enumerate(self.virtual_frames):
            refs = virtual_frame["reflections"]
            refs["intensity.sum.sigma"] = flex.sqrt(refs["intensity.sum.variance"])
            for r in refs.rows():
                h, k, l = r["miller_index"]
                i_sum = r["intensity.sum.value"]
                sig_i_sum = r["intensity.sum.sigma"]
                dfI.loc[frame_id] = [h, k, l, i_sum, sig_i_sum, frame_id + 1]
        dfI["zone_axis_id"] = 0
        lat_params = (a, b, c, alpha, beta, gamma)
        uvws = np.array(uvws)
        self._write_dyn_cif_pets(
            cif_filename, lat_params, volume, wavelength, UBmatrix, uvws, dfI, comment
        )

    def _write_dyn_cif_pets(
        self, name, lat_params, vol, lam, UBmatrix, uvws, dfI, comment=""
    ):
        """
        - name : name of dataset
        - lat_params : (a,b,c,alpha,beta,gamma)
        - vol : volume
        - lam : wavelength
        - UBmatrix : 3x3 matrix orientation
        - uvws : array with columns [u,v,w,scale] (can also contain precession_angle,alpha,beta,omega)
        - dfI : DataFrame intensities ['index_h','index_k','index_l','intensity_meas','intensity_sigma','zone_axis_id']
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
        if not comment:
            comment = """;
        data collection geometry: continuous rotation
        dstarmax:  1.800
        RC width:  0.00120
        mosaicity:  0.210
        rotation axis position:  230.391
        reflection size:    8.000
        Virtual frame settings: number of merged frames:  7
                                step between frames:      5
                                sum all intensities:      1
        ;"""
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
        n_uvws = uvws.shape[0]
        vals = np.hstack([uvws[:, :3], np.zeros((n_uvws, 4)), uvws[:, -1][:, None]])
        for i, v in enumerate(vals):
            lo.add_row(
                [("%d" % (i + 1)).rjust(4)]
                + [("%.5f" % u).rjust(8) for u in v[:3]]
                + ["%.3f" % u for u in v[3:]]
            )

        #### Intensities
        colsI = [
            "index_h",
            "index_k",
            "index_l",
            "intensity_meas",
            "intensity_sigma",
            "zone_axis_id",
        ]
        li = b.init_loop("_refln_", colsI)
        for i, c in dfI.iterrows():
            li.add_row(
                [
                    ("%d" % h).rjust(4)
                    for h in c[["index_h", "index_k", "index_l"]].values
                ]
                + [
                    ("%.2f" % v).rjust(11)
                    for v in c[["intensity_meas", "intensity_sigma"]].values
                ]
                + [("%d" % c["zone_axis_id"]).rjust(5)]
            )

        out = name
        # out = '%s_dyn.cif_files' %name
        doc.write_file(out)
        print("file %s written " % out)
        return out
