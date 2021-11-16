import logging

import numpy as np

from scitbx import matrix

from dials.command_line.frame_orientations import extract_experiment_data

logger = logging.getLogger(__name__)


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
        a table of reflections and the orientation information for the virtual
        frame"""

        # Get the frame orientation data
        directions = self.frame_orientations["directions"]
        zone_axes = self.frame_orientations["zone_axes"]
        real_space_axes = self.frame_orientations["real_space_axes"]
        orientations = self.frame_orientations["orientations"]

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

            # TODO For each reflection, calculate r from s1 and s0. Rotate r according
            # the angle from the centre of the virtual frame. calculate the extinction
            # distance from distance from the centre of the Ewald sphere - |s0|.
            # Actually, according to the Klar et al. paper there are two parameters
            # to consider: Dsg, the distance from the virtual frame edge and Rsg,
            # defined according to the excitation error from the centre of the
            # virtual frame

            # Look up the orientation data using an index, which is the centre
            # of the virtual frame offset so that the scan starts from 0
            centre = int((start + stop) / 2.0)
            index = centre - arr_start

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
