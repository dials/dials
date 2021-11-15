import logging

import numpy as np

logger = logging.getLogger(__name__)


class PETSOutput:
    """Class to export integrated data in PETS CIF format"""

    def __init__(self, experiments, reflections, params):
        """Init with the parameters, experiments and reflections"""

        self.filename_prefix = params.filename_prefix
        self.exp_id = params.id
        self.excitation_error_cutoff = params.virtual_frame.excitation_error_cutoff
        self.n_merged = params.virtual_frame.n_merged
        self.step = params.virtual_frame.step
        self.experiments = experiments
        if len(reflections) != 1:
            raise ValueError("Only a single reflection file can be exported")
        self.reflections = reflections[0]

        self._check_experiments()

    def _check_experiments(self):
        """Check a DIALS experiment list is suitable and convert geometry to the
        PETS coordinate system"""

        if self.exp_id is None:
            if len(self.experiments) != 1:
                raise ValueError("Please select a single experiment for export")
            else:
                self.exp_id = 0

        try:
            experiment = self.experiments[self.exp_id]
        except IndexError as e:
            raise IndexError(
                f"Unable to select experiment {self.exp_id} for export"
            ) from e

        self.reflections = self.reflections.select(
            self.reflections["id"] == self.exp_id
        )

        # Check the unit cell is static
        crystal = experiment.crystal
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
            experiment.goniometer.num_scan_points > 0
            or experiment.beam.num_scan_points > 0
        ):
            raise ValueError("Only scan-varying crystal orientation is supported")

    def write_cif_pets(self):
        pass

    def write_dyn_cif_pets(self):
        pass
