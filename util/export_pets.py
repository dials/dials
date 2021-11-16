import logging
import math
from itertools import tee

import numpy as np

from scitbx import matrix

from dials.array_family import flex

logger = logging.getLogger(__name__)

# An itertools recipe in Python 3.7, but a module function in 3.10
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


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

    def _check_experiments(self):
        """Check a DIALS experiment list is suitable and convert geometry to the
        PETS coordinate system"""

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

        self.reflections = self.reflections.select(
            self.reflections["id"] == self.exp_id
        )

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

        required_keys = ("intensity.sum.value", "intensity.sum.variance")
        if any(x not in self.reflections for x in required_keys):
            raise ValueError(
                f"The reflection table requires these keys: {', '.join(required_keys)}"
            )

        # Select only fully-recorded reflections
        fulls = self.reflections["partiality"] >= self.partiality_cutoff
        logger.info(f"Removing {fulls.count(False)} partial reflections")
        self.reflections = self.reflections.select(fulls)

    def _set_virtual_frames(self):
        """Create a list of virtual frames for the experiment, each containing
        a table of reflections assigned to that virtual frame"""

        # Get reciprocal lattice points from s1 vectors
        s0 = self.experiment.beam.get_s0()
        self.reflections["rlp"] = self.reflections["s1"] - s0

        # Get experimental models
        crystal = self.experiment.crystal
        scan = self.experiment.scan
        goniometer = self.experiment.goniometer
        m2 = matrix.col(goniometer.get_rotation_axis_datum())
        S = matrix.sqr(goniometer.get_setting_rotation())
        F = matrix.sqr(goniometer.get_fixed_rotation())

        # Calculate the orientation matrix in the lab frame for each scan point
        # and set it in the reflections. This is not exactly the U matrix used
        # to integrate the reflections but is close. The individual U matrices
        # for each reflection will be updated afterwards.
        _, _, frames = self.reflections["xyzcal.px"].parts()
        self.reflections["U"] = flex.mat3_double(len(self.reflections))
        Umats = [
            matrix.sqr(crystal.get_U_at_scan_point(i))
            for i in range(crystal.num_scan_points)
        ]
        start, stop = scan.get_array_range()
        SRFUmats = []
        # Loop over the array range frame numbers at the scan-points
        for i, Umat in zip(range(start, stop + 1), Umats):
            phi = scan.get_angle_from_array_index(i, deg=False)
            R = m2.axis_and_angle_as_r3_rotation_matrix(phi, deg=False)
            SRFU = S * R * F * Umat
            SRFUmats.append(SRFU)
            sel = (frames > i) & (frames <= stop)
            self.reflections["U"].set_selected(sel, SRFU)

        # Calculate axis and angle for the transformation at each frame from a
        # comparison of the orientations at the beginning and end of that frame
        ang_ax = []
        for U1, U2 in pairwise(SRFUmats):
            M = U2 * U1.transpose()
            ang_ax.append(
                M.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
                    deg=False
                )
            )

        # Calculate the U matrix for each reflection individually - slow
        for i in range(self.reflections.nrows()):
            ref = self.reflections[i]
            _, _, frame = ref["xyzcal.px"]
            frame_start = math.floor(frame)
            frame_fraction = frame - frame_start
            # Look up the relevant axis and angle for this frame
            lookup = int(frame_start) - start
            angle, axis = ang_ax[lookup]
            # This reflection's position is given by a fractional application of
            # this frame's transformation
            frac_angle = frame_fraction * angle
            Uoffset = axis.axis_and_angle_as_r3_rotation_matrix(frac_angle, deg=False)

            # Update ref["U"] by premultiplying by Uoffset
            self.reflections[i]["U"] = Uoffset * matrix.sqr(ref["U"])

        # Loop over the virtual frames
        starts = list(range(*scan.get_array_range(), self.step))
        for start in starts:
            stop = start + self.n_merged
            centre = (start + stop) / 2.0
            print(centre)

            # Take only reflections whose centroids are inside this virtual frame
            sel = (frames > start) & (frames <= stop)
            refs = self.reflections.select(sel)
            refs["U"] = flex.mat3_double(len(refs))

            pass
        # For each reflection, calculate r from s1 and s0. Rotate r according
        # the angle from the centre of the virtual frame. calculate the extinction
        # distance from distance from the centre of the Ewald sphere - |s0|.

    def write_cif_pets(self):
        pass

    def write_dyn_cif_pets(self):

        self._set_virtual_frames()
